import os
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.linalg import null_space
import matplotlib
from pydantic import BaseModel, Field
from typing import List, Dict, Optional
from fastapi import APIRouter, Response
import redis
import uuid
import io

router = APIRouter()

# ===================== Cylindrical Bending Analysis Code =====================
class Analysis:
    def __init__(self, L, h, q0, angles, material_props, num_layer):
        self.L, self.h, self.q0, self.angles, self.material_props, self.num_layer = L, h, q0, angles, material_props, num_layer
        self.x1, self.x3 = sp.symbols('x1 x3', real=True)
        self.constants = sp.symbols(f'C1:{6*self.num_layer+1}', real=True)

    def solve(self):
        C_matrices = self.compute_C_matrices(self.material_props)
        Rsigma_array = [self.compute_Rsigma_array(self.angles[i]) for i in range(self.num_layer)]
        C_transformed = [self.compute_transformed_C(C_matrices, Rsigma_array[i]) for i in range(self.num_layer)]
        UVWsol = []
        for i in range(self.num_layer):
            GDEsol = self.solve_GDE(C_transformed[i], i)
            UVWsol.append(GDEsol)
        stress0 = [self.get_stress0(C_transformed[i], *UVWsol[i]) for i in range(self.num_layer)]
        eqs = self.determine_BCs_and_continuity(self.q0, UVWsol, stress0, self.num_layer)
        const_sol = self.solve_constants(eqs)
        const_dict = {c: value for c, value in zip(self.constants, const_sol)}
        UVWsol_final = [[expr.subs(const_dict) for expr in sol] for sol in UVWsol]
        displacements = [self.compute_displacement(UVWsol_final[i]) for i in range(self.num_layer)]
        strains = [self.compute_strain(UVWsol_final[i]) for i in range(self.num_layer)]
        stresses = [self.compute_stress(C_transformed[i], strains[i]) for i in range(self.num_layer)]
        displacement_field, strain_field, stress_field = self.compute_layerwise_fields(displacements, strains, stresses)
        return displacement_field, strain_field, stress_field

    def compute_C_matrices(self, material_props):
        E1   = material_props['E1']
        E2   = material_props['E2']
        E3   = material_props['E3']
        G12  = material_props['G12']
        G13  = material_props['G13']
        G23  = material_props['G23']
        nu12 = material_props['nu12']
        nu13 = material_props['nu13']
        nu23 = material_props['nu23']
        Smat = sp.Matrix([
            [1/E1,      -nu12/E1, -nu13/E1, 0,      0,      0],
            [-nu12/E1,  1/E2,     -nu23/E2, 0,      0,      0],
            [-nu13/E1,  -nu23/E2, 1/E3,     0,      0,      0],
            [0,         0,        0,        1/G23,  0,      0],
            [0,         0,        0,        0,      1/G13,  0],
            [0,         0,        0,        0,      0,      1/G12]
        ])
        return Smat.inv()

    def compute_Rsigma_array(self, angle):
        theta = np.pi * angle / 180
        c, s = sp.cos(theta), sp.sin(theta)
        Beta = sp.Matrix([
            [c, -s, 0],
            [s,  c, 0],
            [0,  0, 1]
        ])
        Rsigma = sp.Matrix([
            [Beta[0, 0]**2, Beta[0, 1]**2, Beta[0, 2]**2, 2*Beta[0, 1]*Beta[0, 2], 2*Beta[0, 0]*Beta[0, 2], 2*Beta[0, 0]*Beta[0, 1]],
            [Beta[1, 0]**2, Beta[1, 1]**2, Beta[1, 2]**2, 2*Beta[1, 1]*Beta[1, 2], 2*Beta[1, 0]*Beta[1, 2], 2*Beta[1, 0]*Beta[1, 1]],
            [Beta[2, 0]**2, Beta[2, 1]**2, Beta[2, 2]**2, 2*Beta[2, 1]*Beta[2, 2], 2*Beta[2, 0]*Beta[2, 2], 2*Beta[2, 0]*Beta[2, 1]],
            [Beta[1, 0]*Beta[2, 0], Beta[1, 1]*Beta[2, 1], Beta[1, 2]*Beta[2, 2],
             Beta[1, 2]*Beta[2, 1] + Beta[1, 1]*Beta[2, 2],
             Beta[1, 2]*Beta[2, 0] + Beta[1, 0]*Beta[2, 2],
             Beta[1, 1]*Beta[2, 0] + Beta[1, 0]*Beta[2, 1]],
            [Beta[0, 0]*Beta[2, 0], Beta[0, 1]*Beta[2, 1], Beta[0, 2]*Beta[2, 2],
             Beta[0, 2]*Beta[2, 1] + Beta[0, 1]*Beta[2, 2],
             Beta[0, 2]*Beta[2, 0] + Beta[0, 0]*Beta[2, 2],
             Beta[0, 1]*Beta[2, 0] + Beta[0, 0]*Beta[2, 1]],
            [Beta[0, 0]*Beta[1, 0], Beta[0, 1]*Beta[1, 1], Beta[0, 2]*Beta[1, 2],
             Beta[0, 2]*Beta[1, 1] + Beta[0, 1]*Beta[1, 2],
             Beta[0, 2]*Beta[1, 0] + Beta[0, 0]*Beta[1, 2],
             Beta[0, 1]*Beta[1, 0] + Beta[0, 0]*Beta[1, 1]]
        ])
        return Rsigma

    def compute_transformed_C(self, Cmat, Rsigma):
        return Rsigma * Cmat * Rsigma.T

    def solve_GDE(self, C_transformed, index):
        Lambda = sp.symbols('Lambda')
        Amat = sp.Matrix([
            [-C_transformed[0, 0]*np.pi**2/self.L**2 + C_transformed[4, 4]*Lambda**2,
             -C_transformed[0, 5]*np.pi**2/self.L**2 + C_transformed[3, 4]*Lambda**2,
             C_transformed[0, 2]*np.pi*Lambda/self.L + C_transformed[4, 4]*np.pi*Lambda/self.L],
            [-C_transformed[0, 5]*np.pi**2/self.L**2 + C_transformed[3, 4]*Lambda**2,
             -C_transformed[5, 5]*np.pi**2/self.L**2 + C_transformed[3, 3]*Lambda**2,
             C_transformed[2, 5]*np.pi*Lambda/self.L + C_transformed[3, 4]*np.pi*Lambda/self.L],
            [-C_transformed[0, 2]*np.pi*Lambda/self.L - C_transformed[4, 4]*np.pi*Lambda/self.L,
             -C_transformed[2, 5]*np.pi*Lambda/self.L - C_transformed[3, 4]*np.pi*Lambda/self.L,
             -C_transformed[4, 4]*np.pi**2/self.L**2 + C_transformed[2, 2]*Lambda**2]
        ])
        det_Amat = Amat.det()
        lambda_roots = sp.roots(det_Amat, Lambda)
        lambdasol = [root for root, mult in lambda_roots.items()]
        lambdasol = np.array(lambdasol, dtype=np.complex128)
        null_spaces = []
        for lam in lambdasol:
            Amat_lambda = Amat.subs(Lambda, lam)
            Amat_numeric = np.array(Amat_lambda.evalf(chop=True), dtype=np.complex128)
            null_vectors = null_space(Amat_numeric, rcond=1e-6).flatten()
            null_vectors = np.real_if_close(null_vectors, tol=1e-6)
            null_spaces.append(null_vectors)
        U_expr = sum(self.constants[6*index+i]*null_spaces[i][0]*sp.exp(lambdasol[i]*self.x3) for i in range(6))
        V_expr = sum(self.constants[6*index+i]*null_spaces[i][1]*sp.exp(lambdasol[i]*self.x3) for i in range(6))
        W_expr = sum(self.constants[6*index+i]*null_spaces[i][2]*sp.exp(lambdasol[i]*self.x3) for i in range(6))
        return sp.Matrix([U_expr, V_expr, W_expr])

    def get_stress0(self, C_transformed, U, V, W):
        W_prime = sp.diff(W, self.x3)
        stress0 = sp.zeros(6, 1)
        stress0[0] = -np.pi/self.L * C_transformed[0, 0]*U + C_transformed[0, 2]*W_prime - np.pi/self.L * C_transformed[0, 5]*V
        stress0[1] = -np.pi/self.L * C_transformed[0, 1]*U + C_transformed[1, 2]*W_prime - np.pi/self.L * C_transformed[1, 5]*V
        stress0[2] = -np.pi/self.L * C_transformed[0, 2]*U + C_transformed[2, 2]*W_prime - np.pi/self.L * C_transformed[2, 5]*V
        stress0[3] = C_transformed[3, 3]*sp.diff(V, self.x3) + C_transformed[3, 4]*(sp.diff(U, self.x3)+np.pi/self.L*W)
        stress0[4] = C_transformed[3, 4]*sp.diff(V, self.x3) + C_transformed[4, 4]*(sp.diff(U, self.x3)+np.pi/self.L*W)
        stress0[5] = -np.pi/self.L * C_transformed[0, 5]*U + C_transformed[2, 5]*W_prime - np.pi/self.L * C_transformed[5, 5]*V
        return stress0

    def determine_BCs_and_continuity(self, q0, UVW, stress0, num_layer):
        bcs = []
        continuity = []
        layer_heights = np.linspace(-self.h/2, self.h/2, num_layer+1)
        bcs.append((sp.N(stress0[0][2].subs(self.x3, layer_heights[0])), 0))
        bcs.append((sp.N(stress0[0][3].subs(self.x3, layer_heights[0])), 0))
        bcs.append((sp.N(stress0[0][4].subs(self.x3, layer_heights[0])), 0))
        bcs.append((sp.N(stress0[-1][2].subs(self.x3, layer_heights[-1])), q0))
        bcs.append((sp.N(stress0[-1][3].subs(self.x3, layer_heights[-1])), 0))
        bcs.append((sp.N(stress0[-1][4].subs(self.x3, layer_heights[-1])), 0))
        for i in range(num_layer-1):
            x3_interface = layer_heights[i+1]
            continuity.append((sp.N(UVW[i][0].subs(self.x3, x3_interface) - UVW[i+1][0].subs(self.x3, x3_interface)), 0))
            continuity.append((sp.N(UVW[i][1].subs(self.x3, x3_interface) - UVW[i+1][1].subs(self.x3, x3_interface)), 0))
            continuity.append((sp.N(UVW[i][2].subs(self.x3, x3_interface) - UVW[i+1][2].subs(self.x3, x3_interface)), 0))
            continuity.append((sp.N(stress0[i][2].subs(self.x3, x3_interface) - stress0[i+1][2].subs(self.x3, x3_interface)), 0))
            continuity.append((sp.N(stress0[i][3].subs(self.x3, x3_interface) - stress0[i+1][3].subs(self.x3, x3_interface)), 0))
            continuity.append((sp.N(stress0[i][4].subs(self.x3, x3_interface) - stress0[i+1][4].subs(self.x3, x3_interface)), 0))
        return bcs + continuity

    def solve_constants(self, equations):
        A = []
        b = []
        tol = 1e-6
        for lhs, rhs in equations:
            row = []
            for c in self.constants:
                coeff = lhs.coeff(c)
                coeff = 0 if abs(coeff) < tol else coeff
                row.append(coeff)
            A.append(row)
            b.append(rhs)
        sol = np.linalg.solve(np.array(A, dtype=complex), np.array(b, dtype=complex))
        return sol

    def compute_displacement(self, UVWsol):
        U_i, V_i, W_i = UVWsol
        u1 = U_i * sp.cos(np.pi*self.x1/self.L)
        u2 = V_i * sp.cos(np.pi*self.x1/self.L)
        u3 = W_i * sp.sin(np.pi*self.x1/self.L)
        return sp.Matrix([u1, u2, u3])

    def compute_strain(self, UVWsol):
        U_i, V_i, W_i = UVWsol
        U_prime = sp.diff(U_i, self.x3)
        V_prime = sp.diff(V_i, self.x3)
        W_prime = sp.diff(W_i, self.x3)
        epsilon_11 = -np.pi/self.L * U_i * sp.sin(np.pi*self.x1/self.L)
        epsilon_22 = 0
        epsilon_33 = W_prime * sp.sin(np.pi*self.x1/self.L)
        gamma_23 = V_prime * sp.cos(np.pi*self.x1/self.L)
        gamma_13 = (U_prime + np.pi/self.L*W_i)*sp.cos(np.pi*self.x1/self.L)
        gamma_12 = -np.pi/self.L * V_i * sp.sin(np.pi*self.x1/self.L)
        return sp.Matrix([epsilon_11, epsilon_22, epsilon_33, gamma_23, gamma_13, gamma_12])

    def compute_stress(self, C_transformed, strains):
        return sp.Matrix(C_transformed @ strains)

    def compute_layerwise_fields(self, displacements, strains, stresses):
        layer_boundaries = np.linspace(-self.h/2, self.h/2, self.num_layer+1)
        displacement_piecewise = [[] for _ in range(len(displacements[0]))]
        strain_piecewise = [[] for _ in range(len(strains[0]))]
        stress_piecewise = [[] for _ in range(len(stresses[0]))]
        for i in range(self.num_layer):
            x3_interval = (self.x3 >= layer_boundaries[i]) & (self.x3 < layer_boundaries[i+1])
            for j in range(len(displacements[i])):
                displacement_piecewise[j].append((displacements[i][j], x3_interval))
            for j in range(len(strains[i])):
                strain_piecewise[j].append((strains[i][j], x3_interval))
            for j in range(len(stresses[i])):
                stress_piecewise[j].append((stresses[i][j], x3_interval))
            if i < self.num_layer-1:
                x3_interface = self.x3 == layer_boundaries[i+1]
                for j in range(len(displacements[i])):
                    displacement_piecewise[j].append((displacements[i][j], x3_interface))
                    displacement_piecewise[j].append((displacements[i+1][j], x3_interface))
                for j in range(len(strains[i])):
                    strain_piecewise[j].append((strains[i][j], x3_interface))
                    strain_piecewise[j].append((strains[i+1][j], x3_interface))
                for j in range(len(stresses[i])):
                    stress_piecewise[j].append((stresses[i][j], x3_interface))
                    stress_piecewise[j].append((stresses[i+1][j], x3_interface))
        x3_interval_top = self.x3 == layer_boundaries[-1]
        for j in range(len(displacements[-1])):
            displacement_piecewise[j].append((displacements[-1][j], x3_interval_top))
        for j in range(len(strains[-1])):
            strain_piecewise[j].append((strains[-1][j], x3_interval_top))
        for j in range(len(stresses[-1])):
            stress_piecewise[j].append((stresses[-1][j], x3_interval_top))
        displacement_field = [sp.Piecewise(*dp) for dp in displacement_piecewise]
        strain_field = [sp.Piecewise(*spw) for spw in strain_piecewise]
        stress_field = [sp.Piecewise(*stw) for stw in stress_piecewise]
        return displacement_field, strain_field, stress_field

# ===================== Plotting Function =====================
# Initialize Redis client
def get_redis_client():
    redis_url = os.getenv("REDIS_URL")
    if not redis_url:
        raise ValueError("Failed to get redis url because Heroku Redis is not enabled.")
    return redis.Redis.from_url(redis_url, decode_responses=False)

# Store multiple plots in Redis under a single request ID
def redis_store_plots(plots):
    request_id = str(uuid.uuid4())  # Unique ID per request
    stored_plots = {}

    for plot_name, fig in plots.items():
        buf = io.BytesIO()
        fig.savefig(buf, format="png")  # Save figure as PNG
        buf.seek(0)

        # Store binary data in Redis
        redis_key = f"{request_id}:{plot_name}"
        redis_client = get_redis_client()
        redis_client.setex(redis_key, 3600, buf.getvalue())  # Store for 1 hour (3600s)

        # Store URL reference
        stored_plots[plot_name] = f"/get-plot/{request_id}/{plot_name}"

    return request_id, stored_plots  # Return request ID & URL dictionary

# Use a non-interactive backend for matplotlib
matplotlib.use("Agg")

def plot_results(disp, strain, stress, L, h, x1, x3, disp_x1, strain_x1, stress_x1, plots):
    # Generate 40 equally spaced points along x₁-direction from 0 to L
    x1_vals = np.linspace(0, L, 40)
    
    # Determine the number of sections based on the arguments of the first displacement function
    num_sections = len(disp[0].args)
    base_points = 40 // num_sections
    remainder = 40 % num_sections
    points_per_section = [base_points] * num_sections
    for i in range(remainder):
        points_per_section[i] += 1

    # Define the boundaries of the layers along x₃ from -h/2 to h/2
    layer_boundaries = np.linspace(-h/2, h/2, num_sections+1)
    x3_vals = []
    for i in range(num_sections):
        pts = np.linspace(layer_boundaries[i] + 1e-10, layer_boundaries[i+1] - 1e-10,
                          points_per_section[i], endpoint=True)
        x3_vals.extend(pts)
        if i < num_sections - 1:
            x3_vals.append(layer_boundaries[i+1] - 1e-10)
            x3_vals.append(layer_boundaries[i+1] + 1e-10)
    x3_vals = np.array(x3_vals)
    
    # Prepare grid for full 3D plotting
    X1 = x1_vals[:, None]
    X3 = x3_vals[None, :]

    # Convert symbolic expressions to numerical functions
    disp_funcs = [sp.lambdify((x1, x3), d, "numpy") for d in disp]
    strain_funcs = [sp.lambdify((x1, x3), s, "numpy") for s in strain]
    stress_funcs = [sp.lambdify((x1, x3), s, "numpy") for s in stress]
    
    # Create a dictionary to collect all figures for this request
    all_plots = {}
    
    # If 2D plots are requested
    if "2d_combined" in plots or "2d_standalone" in plots:
        # Evaluate 2D fields at user-specified x₁ values
        disp_numeric_1d = np.array([
            disp_funcs[0](disp_x1[0], x3_vals).real,
            disp_funcs[1](disp_x1[1], x3_vals).real,
            disp_funcs[2](disp_x1[2], x3_vals).real
        ])
        # For strains, here we assume combined order: [ε₁₁, ε₃₃, γ₂₃, γ₁₃, γ₁₂]
        strain_numeric_1d = np.array([
            strain_funcs[0](strain_x1[0], x3_vals).real,
            strain_funcs[2](strain_x1[1], x3_vals).real,
            strain_funcs[3](strain_x1[2], x3_vals).real,
            strain_funcs[4](strain_x1[3], x3_vals).real,
            strain_funcs[5](strain_x1[4], x3_vals).real
        ])
        stress_numeric_1d = np.array([
            stress_funcs[0](stress_x1[0], x3_vals).real,
            stress_funcs[1](stress_x1[1], x3_vals).real,
            stress_funcs[2](stress_x1[2], x3_vals).real,
            stress_funcs[3](stress_x1[3], x3_vals).real,
            stress_funcs[4](stress_x1[4], x3_vals).real,
            stress_funcs[5](stress_x1[5], x3_vals).real
        ])
        disp_labels = [
            rf"$U$ at $x_1 = {disp_x1[0]} (m)$",
            rf"$V$ at $x_1 = {disp_x1[1]} (m)$",
            rf"$W$ at $x_1 = {disp_x1[2]} (m)$"
        ]
        strain_labels = [
            rf"$\epsilon_{{11}}$ at $x_1 = {strain_x1[0]} (m)$",
            rf"$\epsilon_{{33}}$ at $x_1 = {strain_x1[1]} (m)$",
            rf"$\gamma_{{23}}$ at $x_1 = {strain_x1[2]} (m)$",
            rf"$\gamma_{{13}}$ at $x_1 = {strain_x1[3]} (m)$",
            rf"$\gamma_{{12}}$ at $x_1 = {strain_x1[4]} (m)$"
        ]
        stress_labels = [
            rf"$\sigma_{{11}}$ at $x_1 = {stress_x1[0]} (m)$",
            rf"$\sigma_{{22}}$ at $x_1 = {stress_x1[1]} (m)$",
            rf"$\sigma_{{33}}$ at $x_1 = {stress_x1[2]} (m)$",
            rf"$\sigma_{{23}}$ at $x_1 = {stress_x1[3]} (m)$",
            rf"$\sigma_{{13}}$ at $x_1 = {stress_x1[4]} (m)$",
            rf"$\sigma_{{12}}$ at $x_1 = {stress_x1[5]} (m)$"
        ]

        # If 2D combined plots are requested
        if "2d_combined" in plots:
            # Combined 2D Displacement Plot
            fig_disp = plt.figure(figsize=(6, 4))
            for i in range(len(disp_numeric_1d)):
                plt.plot(x3_vals, disp_numeric_1d[i], linewidth=2, label=disp_labels[i])
            plt.xlabel(r"$x_3$ (m)")
            plt.ylabel("Displacement (m)")
            plt.legend()
            plt.grid(True)
            fig_disp.tight_layout()
            all_plots["2d_displacement"] = fig_disp
            plt.close(fig_disp)

            # Combined 2D Strain Plot
            fig_strain = plt.figure(figsize=(6, 4))
            for i in range(len(strain_numeric_1d)):
                plt.plot(x3_vals, strain_numeric_1d[i], linewidth=2, label=strain_labels[i])
            plt.xlabel(r"$x_3$ (m)")
            plt.ylabel("Strain")
            plt.legend()
            plt.grid(True)
            fig_strain.tight_layout()
            all_plots["2d_strain"] = fig_strain
            plt.close(fig_strain)

            # Combined 2D Stress Plot
            fig_stress = plt.figure(figsize=(6, 4))
            for i in range(len(stress_numeric_1d)):
                plt.plot(x3_vals, stress_numeric_1d[i], linewidth=2, label=stress_labels[i])
            plt.xlabel(r"$x_3$ (m)")
            plt.ylabel("Stress (Pa)")
            plt.legend()
            plt.grid(True)
            fig_stress.tight_layout()
            all_plots["2d_stress"] = fig_stress
            plt.close(fig_stress)

        # If 2D separate plots are requested
        if "2d_standalone" in plots:
            # Collect standalone 2D Displacement Plots
            for i in range(len(disp_numeric_1d)):
                fig = plt.figure(figsize=(6, 4))
                plt.plot(x3_vals, disp_numeric_1d[i], linewidth=2, label=disp_labels[i])
                plt.xlabel(r"$x_3$ (m)")
                plt.ylabel("Displacement (m)")
                plt.legend()
                plt.grid(True)
                all_plots[f"2d_disp_standalone_{i}"] = fig
                plt.close(fig)

            # Collect standalone 2D Strain Plots
            for i in range(len(strain_numeric_1d)):
                fig = plt.figure(figsize=(6, 4))
                plt.plot(x3_vals, strain_numeric_1d[i], linewidth=2, label=strain_labels[i])
                plt.xlabel(r"$x_3$ (m)")
                plt.ylabel("Strain")
                plt.legend()
                plt.grid(True)
                all_plots[f"2d_strain_standalone_{i}"] = fig
                plt.close(fig)

            # Collect standalone 2D Stress Plots
            for i in range(len(stress_numeric_1d)):
                fig = plt.figure(figsize=(6, 4))
                plt.plot(x3_vals, stress_numeric_1d[i], linewidth=2, label=stress_labels[i])
                plt.xlabel(r"$x_3$ (m)")
                plt.ylabel("Stress (Pa)")
                plt.legend()
                plt.grid(True)
                all_plots[f"2d_stress_standalone_{i}"] = fig
                plt.close(fig)

    # If 3D plots are requested
    if "3d" in plots:
        # Prepare grid for 3D plots
        def evaluate_functions(func_list, threshold):
            results = np.array([func(X1, X3).real for func in func_list], dtype=object)
            for i in range(len(results)):
                results[i][np.abs(results[i]) < threshold] = 0
            return results

        disp_numeric_3d = evaluate_functions(disp_funcs, 1e-8)
        strain_numeric_3d = evaluate_functions(strain_funcs, 1e-10)
        stress_numeric_3d = evaluate_functions(stress_funcs, 1e-4)

        # 3D Displacement Plots
        disp_titles_3d = [r"$u_1$ (m)", r"$u_2$ (m)", r"$u_3$ (m)"]
        for i in range(len(disp_numeric_3d)):
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X1, X3, disp_numeric_3d[i], cmap='gist_rainbow')
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_3$")
            ax.set_zlabel(disp_titles_3d[i])
            ax.set_title(f"Displacement {disp_titles_3d[i]}")
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
            fig.tight_layout()
            all_plots[f"3d_displacement_{i}"] = fig  # Collect figure
            plt.close(fig)

        # 3D Strain Plots
        strain_titles_3d = [r"$\epsilon_{11}$", r"$\epsilon_{22}$", r"$\epsilon_{33}$",
                            r"$\gamma_{23}$", r"$\gamma_{13}$", r"$\gamma_{12}$"]
        for i in range(len(strain_numeric_3d)):
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X1, X3, strain_numeric_3d[i], cmap='gist_rainbow')
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_3$")
            ax.set_zlabel(strain_titles_3d[i])
            ax.set_title(f"Strain {strain_titles_3d[i]}")
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
            fig.tight_layout()
            all_plots[f"3d_strain_{i}"] = fig
            plt.close(fig)

        # 3D Stress Plots
        stress_titles_3d = [r"$\sigma_{11}$ (Pa)", r"$\sigma_{22}$ (Pa)", r"$\sigma_{33}$ (Pa)",
                            r"$\sigma_{23}$ (Pa)", r"$\sigma_{13}$ (Pa)", r"$\sigma_{12}$ (Pa)"]
        for i in range(len(stress_numeric_3d)):
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X1, X3, stress_numeric_3d[i], cmap='gist_rainbow')
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_3$")
            ax.set_zlabel(stress_titles_3d[i])
            ax.set_title(f"Stress {stress_titles_3d[i]}")
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
            fig.tight_layout()
            all_plots[f"3d_stress_{i}"] = fig
            plt.close(fig)
    
    # Store in redis after collecting all requested plots
    plots_id, figures = redis_store_plots(all_plots)  # Store plots and get URLs

    # Return the request_id and figures (URLs) in the response
    return {"plots_id": plots_id, "figures": figures}

# ===================== Probing Function =====================
def probe_value(disp, strain, stress, L, h, x1_input, x3_input, x1, x3):
    try:
        if not (0 <= x1_input <= L):
            return {"error": f"x1={x1_input} is out of bounds. It should be within [0, {L}]."}
        if not (-h/2 <= x3_input <= h/2):
            return {"error": f"x3={x3_input} is out of bounds. It should be within [{-h/2}, {h/2}]."}
        
        disp_funcs = [sp.lambdify((x1, x3), d, "numpy") for d in disp]
        strain_funcs = [sp.lambdify((x1, x3), s, "numpy") for s in strain]
        stress_funcs = [sp.lambdify((x1, x3), s, "numpy") for s in stress]

        disp_probe = {
            "U": float(disp_funcs[0](x1_input, x3_input).real),
            "V": float(disp_funcs[1](x1_input, x3_input).real),
            "W": float(disp_funcs[2](x1_input, x3_input).real)
        }
        strain_probe = {
            "ε_11": float(strain_funcs[0](x1_input, x3_input).real),
            "ε_22": float(strain_funcs[1](x1_input, x3_input).real),
            "ε_33": float(strain_funcs[2](x1_input, x3_input).real),
            "γ_23": float(strain_funcs[3](x1_input, x3_input).real),
            "γ_13": float(strain_funcs[4](x1_input, x3_input).real),
            "γ_12": float(strain_funcs[5](x1_input, x3_input).real)
        }
        stress_probe = {
            "σ_11": float(stress_funcs[0](x1_input, x3_input).real),
            "σ_22": float(stress_funcs[1](x1_input, x3_input).real),
            "σ_33": float(stress_funcs[2](x1_input, x3_input).real),
            "σ_23": float(stress_funcs[3](x1_input, x3_input).real),
            "σ_13": float(stress_funcs[4](x1_input, x3_input).real),
            "σ_12": float(stress_funcs[5](x1_input, x3_input).real)
        }
        return {"displacement": disp_probe, "strain": strain_probe, "stress": stress_probe}

    except ValueError as e:
        return {"error": str(e)}
    except Exception as e:
        return {"error": f"Error in probing: {str(e)}"}

# ===================== Run analysis =====================
def run_laminate_analysis(L, h, q0, layer_angles, material_props, disp_x1=None, strain_x1=None, stress_x1=None, probe_x1=None, probe_x3=None, plots=None):
    input_error = check_inputs(L, h, layer_angles, material_props)
    if input_error:
        return {"error": input_error}
    if disp_x1 is None:
        disp_x1 = [0, 0, L/2]
    if strain_x1 is None:
        strain_x1 = [L/2, L/2, 0, 0, L/2]
    if stress_x1 is None:
        stress_x1 = [L/2, L/2, L/2, 0, 0, L/2]
    # Perform analysis
    analysis = Analysis(L, h, q0, layer_angles, material_props, len(layer_angles))
    # Solve for displacement, strain, and stress
    disp, strain, stress = analysis.solve()
    # Initialize results as a dictionary
    result_dict = {}
    # Generate plots and get the id assigned to the plots
    plots_dict = plot_results(disp, strain, stress, L, h, analysis.x1, analysis.x3, disp_x1, strain_x1, stress_x1, plots)
    if plots_dict:
        result_dict["plots_id"] = plots_dict["plots_id"]
        result_dict["figures"] = plots_dict["figures"]
    # Probe values if x1 and x3 are specified
    if probe_x1 is not None and probe_x3 is not None:
        probed_results = probe_value(disp, strain, stress, L, h, probe_x1, probe_x3, analysis.x1, analysis.x3)
        result_dict["probe_results"] = probed_results

    return result_dict

# ===================== Check inputs =====================
def check_inputs(L, h, layer_angles, material_props):
    # Check geometry
    if not L > 0:
        return f"Length {L} must be positive."
    if not h > 0:
        return f"Height {h} must be positive."
    if not len(layer_angles) > 0:
        return f"The number of layers {len(layer_angles)} must be a positive integer."
    
    # Ensure all required properties exist
    required_keys = ['E1', 'E2', 'E3', 'G12', 'G13', 'G23', 'nu12', 'nu13', 'nu23']
    missing_keys = [key for key in required_keys if key not in material_props]
    if missing_keys:
        return f"Missing material properties: {', '.join(missing_keys)}."

    # Check positive values for stiffness parameters
    invalid_moduli = [k for k in ['E1', 'E2', 'E3', 'G12', 'G13', 'G23'] if material_props[k] <= 0]
    if invalid_moduli:
        return f"The following parameters must be positive: {', '.join(invalid_moduli)}."

    # Assign values
    E1, E2, E3 = material_props['E1'], material_props['E2'], material_props['E3']
    nu12, nu13, nu23 = material_props['nu12'], material_props['nu13'], material_props['nu23']
    # Compute Poisson's ratios
    nu21 = (E2 / E1) * nu12
    nu31 = (E3 / E1) * nu13
    nu32 = (E3 / E2) * nu23
    # Check material stability conditions
    if not (1 - nu12 * nu21 > 0 and
            1 - nu13 * nu31 > 0 and
            1 - nu23 * nu32 > 0 and
            1 - nu12 * nu21 - nu23 * nu32 - nu13 * nu31 - 2 * nu21 * nu13 * nu32 > 0):
        return "Given material properties are unrealistic."

    return None  # Return None if all inputs are valid

# ===================== Cylindrical Bending Analysis with Input/Output =====================
class CylindricalBendingInput(BaseModel):
    L: float = Field(4.0, description="Length of the laminate (m)")
    h: float = Field(1.0, description="Thickness of the laminate (m)")
    q0: float = Field(1e6, description="Applied load on the laminate (Pa)")
    layer_angles: List[float] = Field([0], description="List of layer angles in degrees")
    material_props: Dict[str, float] = Field({
        "E1": 140000000000, "E2": 10000000000, "E3": 10000000000,
        "G12": 7000000000, "G13": 7000000000, "G23": 3360000000,
        "nu12": 0.3, "nu13": 0.3, "nu23": 0.49
    }, description="Material properties")
    disp_x1: Optional[List[float]] = Field(None, description="x1 coordinates for displacement evaluation")
    strain_x1: Optional[List[float]] = Field(None, description="x1 coordinates for strain evaluation")
    stress_x1: Optional[List[float]] = Field(None, description="x1 coordinates for stress evaluation")
    probe_x1: Optional[float] = Field(None, description="x1 coordinate for probing specific values")
    probe_x3: Optional[float] = Field(None, description="x3 coordinate for probing specific values")
    plots: List[str] = Field(["2d_combined"], description="List of plot types to generate (e.g., '2d_combined', '3d')")

class CylindricalBendingOutput(BaseModel):
    figures: Dict[str, str] = Field(default_factory=dict, description="Mapping of plot names to image URLs")
    probe_results: Optional[Dict] = Field(None, description="Probed displacement, strain, and stress values if provided")
    error: Optional[str] = Field(None, description="Error message if an error occurred")

# ----------------- Handle API calls -----------------
@router.post("/cylindrical-bending-analysis", response_model=CylindricalBendingOutput)
def laminate_analysis(request: CylindricalBendingInput):
    """
    Perform laminate analysis based on provided parameters.

    - **L**: Length of the laminate.
    - **h**: Thickness of the laminate.
    - **q0**: Applied load.
    - **layer_angles**: List of angles for each layer (in degrees).
    - **material_props**: Dictionary containing material properties.
    - **disp_x1**, **strain_x1**, **stress_x1**: x1 coordinates for evaluating fields.
    - **probe_x1**, **probe_x3**: Coordinates for probing specific values.
    - **plots**: List specifying which plots to generate.

    Returns a JSON object containing URLs for generated figures and, if applicable, the probed numerical results.

    ### Example `curl` Request:
    ```bash
    curl -X 'POST' \\
      'http://127.0.0.1:8000/cylindrical-bending-analysis' \\
      -H 'accept: application/json' \\
      -H 'Content-Type: application/json' \\
      -d '{
        "L": 4.0,
        "h": 1.0,
        "q0": 1000000,
        "layer_angles": [0, 45],
        "material_props": {
          "E1": 140000000000,
          "E2": 10000000000,
          "E3": 10000000000,
          "G12": 7000000000,
          "G13": 7000000000,
          "G23": 3360000000,
          "nu12": 0.3,
          "nu13": 0.3,
          "nu23": 0.49
        },
        "disp_x1": [0, 0, 2],
        "strain_x1": [2, 2, 0, 0, 2],
        "stress_x1": [2, 2, 2, 0, 0, 2],
        "plots": ["2d_combined", "3d"]
      }'
    ```
    """
    try:
        results_dict = run_laminate_analysis(
            L=request.L, h=request.h, q0=request.q0,
            layer_angles=request.layer_angles,
            material_props=request.material_props,
            disp_x1=request.disp_x1,
            strain_x1=request.strain_x1,
            stress_x1=request.stress_x1,
            probe_x1=request.probe_x1,
            probe_x3=request.probe_x3,
            plots = request.plots
        )
        # If an error occurred during analysis, return it in the output.
        if "error" in results_dict:
            return CylindricalBendingOutput(figures={}, probe_results=None, error=results_dict["error"])
    except ValueError as e:
        return CylindricalBendingOutput(figures={}, probe_results=None, error=str(e))
    except Exception as e:
        return CylindricalBendingOutput(figures={}, probe_results=None, error=f"Unexpected error: {str(e)}")
    
    # Get figures from results_dict
    figures = results_dict.get("figures", {})

    # Get probe results
    probe_results = results_dict.get("probe_results", None)

    return CylindricalBendingOutput(figures=figures, probe_results=probe_results, error = None)

# Retrieve multiple plots from Redis using the request ID
@router.get("/get-plot/{request_id}/{plot_name}")
def get_plot(request_id: str, plot_name: str):
    redis_key = f"{request_id}:{plot_name}"
    redis_client = get_redis_client()
    image_data = redis_client.get(redis_key)

    if image_data is None:
        return Response(content="Plot not found or expired.", status_code=404)

    return Response(content=image_data, media_type="image/png")
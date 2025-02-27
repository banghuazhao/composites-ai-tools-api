import os
import json
import time
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.linalg import null_space
import matplotlib
from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Union
from fastapi import APIRouter, HTTPException

router = APIRouter()

# Use a non-interactive backend for matplotlib
matplotlib.use("Agg")

# Ensure "results" directory exists
RESULTS_DIR = os.path.join(os.getcwd(), "results")
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# ----------------- Helper Function to Get ngrok URL -----------------
def get_backend_url():
    return os.environ.get("BACKEND_URL", "http://127.0.0.1:8000")

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
def plot_results(disp, strain, stress, L, h, x1, x3, disp_x1, strain_x1, stress_x1, plots):
    # --- 2D Plots Setup (common to both combined and standalone) ---
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
                          points_per_section[i], endpoint=False)
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

        # --- Run only requested plots ---
        if "2d_combined" in plots:
            # Combined 2D Displacement Plot
            plt.figure(figsize=(8, 6))
            for i in range(len(disp_numeric_1d)):
                plt.plot(x3_vals, disp_numeric_1d[i], linewidth=2, label=disp_labels[i])
            plt.xlabel(r"$x_3$ (m)")
            plt.ylabel("Displacement (m)")
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(RESULTS_DIR, "fig2d_displacement.png"))
            plt.close()

            # Combined 2D Strain Plot
            plt.figure(figsize=(8, 6))
            for i in range(len(strain_numeric_1d)):
                plt.plot(x3_vals, strain_numeric_1d[i], linewidth=2, label=strain_labels[i])
            plt.xlabel(r"$x_3$ (m)")
            plt.ylabel("Strain")
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(RESULTS_DIR, "fig2d_strain.png"))
            plt.close()

            # Combined 2D Stress Plot
            plt.figure(figsize=(8, 6))
            for i in range(len(stress_numeric_1d)):
                plt.plot(x3_vals, stress_numeric_1d[i], linewidth=2, label=stress_labels[i])
            plt.xlabel(r"$x_3$ (m)")
            plt.ylabel("Stress (Pa)")
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(RESULTS_DIR, "fig2d_stress.png"))
            plt.close()

        if "2d_standalone" in plots:
            # Standalone 2D Displacement Plots
            for i in range(len(disp_numeric_1d)):
                plt.figure(figsize=(8, 6))
                plt.plot(x3_vals, disp_numeric_1d[i], linewidth=2, label=disp_labels[i])
                plt.xlabel(r"$x_3$ (m)")
                plt.ylabel("Displacement (m)")
                plt.legend()
                plt.grid(True)
                plt.savefig(os.path.join(RESULTS_DIR, f"fig2d_disp_{i}.png"))
                plt.close()

            # Standalone 2D Strain Plots
            for i in range(len(strain_numeric_1d)):
                plt.figure(figsize=(8, 6))
                plt.plot(x3_vals, strain_numeric_1d[i], linewidth=2, label=strain_labels[i])
                plt.xlabel(r"$x_3$ (m)")
                plt.ylabel("Strain")
                plt.legend()
                plt.grid(True)
                plt.savefig(os.path.join(RESULTS_DIR, f"fig2d_strain_{i}.png"))
                plt.close()

            # Standalone 2D Stress Plots
            for i in range(len(stress_numeric_1d)):
                plt.figure(figsize=(8, 6))
                plt.plot(x3_vals, stress_numeric_1d[i], linewidth=2, label=stress_labels[i])
                plt.xlabel(r"$x_3$ (m)")
                plt.ylabel("Stress (Pa)")
                plt.legend()
                plt.grid(True)
                plt.savefig(os.path.join(RESULTS_DIR, f"fig2d_stress_{i}.png"))
                plt.close()

    if "3d" in plots:
        # Prepare grid for 3D plots (reuse X1, X3)
        def evaluate_functions(func_list, threshold):
            results = np.array([func(X1, X3).real for func in func_list], dtype=object)
            for i in range(len(results)):
                results[i][np.abs(results[i]) < threshold] = 0
            return results

        disp_numeric_3d = evaluate_functions(disp_funcs, 1e-8)
        strain_numeric_3d = evaluate_functions(strain_funcs, 1e-10)
        stress_numeric_3d = evaluate_functions(stress_funcs, 1e-4)

        # 3D Displacement Plots (Standalone per component)
        disp_titles_3d = [r"$u_1$ (m)", r"$u_2$ (m)", r"$u_3$ (m)"]
        for i in range(len(disp_numeric_3d)):
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X1, X3, disp_numeric_3d[i], cmap='gist_rainbow')
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_3$")
            ax.set_zlabel(disp_titles_3d[i])
            ax.set_title(f"3D Displacement {disp_titles_3d[i]}")
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
            fig.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, f"fig3d_disp_{i}.png"))
            plt.close(fig)

        # 3D Strain Plots (Standalone per component)
        strain_titles_3d = [r"$\epsilon_{11}$", r"$\epsilon_{22}$", r"$\epsilon_{33}$",
                            r"$\gamma_{23}$", r"$\gamma_{13}$", r"$\gamma_{12}$"]
        for i in range(len(strain_numeric_3d)):
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X1, X3, strain_numeric_3d[i], cmap='gist_rainbow')
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_3$")
            ax.set_zlabel(strain_titles_3d[i])
            ax.set_title(f"3D Strain {strain_titles_3d[i]}")
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
            fig.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, f"fig3d_strain_{i}.png"))
            plt.close(fig)

        # 3D Stress Plots (Standalone per component)
        stress_titles_3d = [r"$\sigma_{11}$ (Pa)", r"$\sigma_{22}$ (Pa)", r"$\sigma_{33}$ (Pa)",
                            r"$\sigma_{23}$ (Pa)", r"$\sigma_{13}$ (Pa)", r"$\sigma_{12}$ (Pa)"]
        for i in range(len(stress_numeric_3d)):
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X1, X3, stress_numeric_3d[i], cmap='gist_rainbow')
            ax.set_xlabel("$x_1$")
            ax.set_ylabel("$x_3$")
            ax.set_zlabel(stress_titles_3d[i])
            ax.set_title(f"3D Stress {stress_titles_3d[i]}")
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
            fig.tight_layout()
            plt.savefig(os.path.join(RESULTS_DIR, f"fig3d_stress_{i}.png"))
            plt.close(fig)

# ===================== Probing Functions =====================
def probe_value(disp, strain, stress, L, h, x1_input, x3_input, x1, x3):
    try:
        if not (0 <= x1_input <= L):
            raise ValueError(f"x1={x1_input} is out of bounds. It should be within [0, {L}].")
        if not (-h/2 <= x3_input <= h/2):
            raise ValueError(f"x3={x3_input} is out of bounds. It should be within [{-h/2}, {h/2}].")
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
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error in probing: {str(e)}")

# ===================== Configuration =====================
def load_config(config_file="config.json"):
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            return json.load(f)
    return {}

def run_laminate_analysis(L=4.0, h=1.0, q0=1e6, layer_angles=[0],
                          material_props=None, disp_x1=None, strain_x1=None,
                          stress_x1=None, probe_x1=None, probe_x3=None, plots=["3d"]):
    config = load_config("config.json")
    L = config.get("L", L)
    h = config.get("h", h)
    q0 = config.get("q0", q0)
    layer_angles = config.get("layer_angles", layer_angles)
    if material_props is None:
        material_props = config.get("material_props", {
            'E1': 140e9, 'E2': 10e9, 'E3': 10e9,
            'G12': 7e9, 'G13': 7e9, 'G23': 3.36e9,
            'nu12': 0.3, 'nu13': 0.3, 'nu23': 0.49
        })
    if disp_x1 is None:
        disp_x1 = config.get("disp_x1", [0, 0, L/2])
    if strain_x1 is None:
        strain_x1 = config.get("strain_x1", [L/2, L/2, 0, 0, L/2])
    if stress_x1 is None:
        stress_x1 = config.get("stress_x1", [L/2, L/2, L/2, 0, 0, L/2])
    num_layers = len(layer_angles)
    analysis = Analysis(L, h, q0, layer_angles, material_props, num_layers)
    disp, strain, stress = analysis.solve()
    plot_results(disp, strain, stress, L, h, analysis.x1, analysis.x3, disp_x1, strain_x1, stress_x1, plots)
    if probe_x1 is not None and probe_x3 is not None:
        probed_results = probe_value(disp, strain, stress, L, h, probe_x1, probe_x3, analysis.x1, analysis.x3)
        return {"probe_results": probed_results}
    return {}

# ===================== Cylindrical Bending Analysis with Input/Output =====================
class CylindricalBendingInput(BaseModel):
    L: float = Field(4.0, description="Length of the laminate (m)")
    h: float = Field(1.0, description="Thickness of the laminate (m)")
    q0: float = Field(1e6, description="Applied load on the laminate (Pa)")
    layer_angles: List[float] = Field([0], description="List of layer angles in degrees")
    material_props: Optional[Dict] = Field(None, description="Material properties (e.g., E1, E2, E3, G12, G13, G23, nu12, nu13, nu23)")
    disp_x1: Optional[List[float]] = Field(None, description="x1 coordinates for displacement evaluation")
    strain_x1: Optional[List[float]] = Field(None, description="x1 coordinates for strain evaluation")
    stress_x1: Optional[List[float]] = Field(None, description="x1 coordinates for stress evaluation")
    probe_x1: Optional[float] = Field(None, description="x1 coordinate for probing specific values")
    probe_x3: Optional[float] = Field(None, description="x3 coordinate for probing specific values")
    plots: Optional[List[str]] = Field(None, description="List of plot types to generate (e.g., '2d_combined', '3d')")

class CylindricalBendingOutput(BaseModel):
    figures: Dict[str, Union[str, List[str]]] = Field(default_factory=dict, description="Mapping of plot names to accessible URLs or list of URLs")
    probe_results: Optional[Dict] = Field(None, description="Probed displacement, strain, and stress values if provided")

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
        results = run_laminate_analysis(
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
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
    
    public_url = get_backend_url()
    timestamp = int(time.time())  # Generate a unique timestamp
    figures = {}
    probe_results = None

    if request.plots and "2d_combined" in request.plots:
        figures.update({
            "2d_displacement": f"{public_url}/results/fig2d-displacement.png",
            "2d_strain": f"{public_url}/results/fig2d-strain.png",
            "2d_stress": f"{public_url}/results/fig2d-stress.png"
        })

    if request.plots and "2d_standalone" in request.plots:
        figures.update({
            "2d_displacement": [f"{public_url}/results/fig2d-disp-{i}.png" for i in range(3)],
            "2d_strain": [f"{public_url}/results/fig2d-strain-{i}.png" for i in range(5)],
            "2d_stress": [f"{public_url}/results/fig2d-stress-{i}.png" for i in range(6)]
        })

    if request.plots and "3d" in request.plots:
        figures.update({
            "3d_displacement": [f"{public_url}/results/fig3d-disp-{i}.png" for i in range(3)],
            "3d_strain": [f"{public_url}/results/fig3d-strain-{i}.png" for i in range(6)],
            "3d_stress": [f"{public_url}/results/fig3d-stress-{i}.png" for i in range(6)]
        })

    if request.plots and "probe" in request.plots and isinstance(results, dict) and "probe_results" in results:
        probe_results = {
            "probe_displacement": results["probe_results"]["displacement"],
            "probe_strain": results["probe_results"]["strain"],
            "probe_stress": results["probe_results"]["stress"]
        }

    return CylindricalBendingOutput(figures=figures, probe_results=probe_results)
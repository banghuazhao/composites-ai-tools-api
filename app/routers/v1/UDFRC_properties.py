from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, Tuple
import numpy as np
from app.models.v1.three_dimensional_properties_output import ThreeDimensionalPropertiesOutput

router = APIRouter()

class FiberProperties(BaseModel):
    E1: float
    E2: float
    G12: float
    nu12: float
    nu23: float

class MatrixProperties(BaseModel):
    E1: float
    nu: float

class UDFRCProperties(BaseModel):
    fiber_properties: FiberProperties
    matrix_properties: MatrixProperties
    fiber_volume_fraction: float

def calculate_voigt_properties(fiber_props: FiberProperties, matrix_props: MatrixProperties, fiber_volume_fraction: float) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    Vf = fiber_volume_fraction
    Vm = 1 - Vf

    ef1 = fiber_props.E1
    ef2 = fiber_props.E2
    ef3 = fiber_props.E2
    gf12 = fiber_props.G12
    gf13 = fiber_props.G12
    gf23 = fiber_props.E2 / (2 * (1 + fiber_props.nu23))
    nuf12 = fiber_props.nu12
    nuf13 = fiber_props.nu12
    nuf23 = fiber_props.nu23

    em1 = matrix_props.E1
    em2 = matrix_props.E1
    em3 = matrix_props.E1
    gm12 = matrix_props.E1 / (2 * (1 + matrix_props.nu))
    gm13 = gm12
    gm23 = gm12
    num12 = matrix_props.nu
    num13 = matrix_props.nu
    num23 = matrix_props.nu

    Sf = np.array([
        [1 / ef1, -nuf12 / ef1, -nuf13 / ef1, 0, 0, 0],
        [-nuf12 / ef1, 1 / ef2, -nuf23 / ef2, 0, 0, 0],
        [-nuf13 / ef1, -nuf23 / ef2, 1 / ef3, 0, 0, 0],
        [0, 0, 0, 1 / gf23, 0, 0],
        [0, 0, 0, 0, 1 / gf13, 0],
        [0, 0, 0, 0, 0, 1 / gf12]
    ])
    
    Sm = np.array([
        [1 / em1, -num12 / em1, -num13 / em1, 0, 0, 0],
        [-num12 / em1, 1 / em2, -num23 / em2, 0, 0, 0],
        [-num13 / em1, -num23 / em2, 1 / em3, 0, 0, 0],
        [0, 0, 0, 1 / gm23, 0, 0],
        [0, 0, 0, 0, 1 / gm13, 0],
        [0, 0, 0, 0, 0, 1 / gm12]
    ])
    
    Cf = np.linalg.inv(Sf)
    Cm = np.linalg.inv(Sm)
    
    CVs = Cf * Vf + Cm * Vm
    SVs = np.linalg.inv(CVs)

    engineering_constants_voigt = {
        "E1": 1 / SVs[0, 0],
        "E2": 1 / SVs[1, 1],
        "E3": 1 / SVs[2, 2],
        "G12": 1 / SVs[5, 5],
        "G13": 1 / SVs[4, 4],
        "G23": 1 / SVs[3, 3],
        "nu12": -SVs[0, 1] / SVs[0, 0],
        "nu13": -SVs[0, 2] / SVs[0, 0],
        "nu23": -SVs[1, 2] / SVs[1, 1]
    }

    return CVs, SVs, engineering_constants_voigt

def calculate_reuss_properties(fiber_props: FiberProperties, matrix_props: MatrixProperties, fiber_volume_fraction: float) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    Vf = fiber_volume_fraction
    Vm = 1 - Vf

    ef1 = fiber_props.E1
    ef2 = fiber_props.E2
    ef3 = fiber_props.E2
    gf12 = fiber_props.G12
    gf13 = fiber_props.G12
    gf23 = fiber_props.E2 / (2 * (1 + fiber_props.nu23))
    nuf12 = fiber_props.nu12
    nuf13 = fiber_props.nu12
    nuf23 = fiber_props.nu23

    em1 = matrix_props.E1
    em2 = matrix_props.E1
    em3 = matrix_props.E1
    gm12 = matrix_props.E1 / (2 * (1 + matrix_props.nu))
    gm13 = gm12
    gm23 = gm12
    num12 = matrix_props.nu
    num13 = matrix_props.nu
    num23 = matrix_props.nu

    Sf = np.array([
        [1 / ef1, -nuf12 / ef1, -nuf13 / ef1, 0, 0, 0],
        [-nuf12 / ef1, 1 / ef2, -nuf23 / ef2, 0, 0, 0],
        [-nuf13 / ef1, -nuf23 / ef2, 1 / ef3, 0, 0, 0],
        [0, 0, 0, 1 / gf23, 0, 0],
        [0, 0, 0, 0, 1 / gf13, 0],
        [0, 0, 0, 0, 0, 1 / gf12]
    ])
    
    Sm = np.array([
        [1 / em1, -num12 / em1, -num13 / em1, 0, 0, 0],
        [-num12 / em1, 1 / em2, -num23 / em2, 0, 0, 0],
        [-num13 / em1, -num23 / em2, 1 / em3, 0, 0, 0],
        [0, 0, 0, 1 / gm23, 0, 0],
        [0, 0, 0, 0, 1 / gm13, 0],
        [0, 0, 0, 0, 0, 1 / gm12]
    ])
    
    SRs = Sf * Vf + Sm * Vm
    CRs = np.linalg.inv(SRs)

    engineering_constants_reuss = {
        "E1": 1 / SRs[0, 0],
        "E2": 1 / SRs[1, 1],
        "E3": 1 / SRs[2, 2],
        "G12": 1 / SRs[5, 5],
        "G13": 1 / SRs[4, 4],
        "G23": 1 / SRs[3, 3],
        "nu12": -SRs[0, 1] / SRs[0, 0],
        "nu13": -SRs[0, 2] / SRs[0, 0],
        "nu23": -SRs[1, 2] / SRs[1, 1]
    }

    return CRs, SRs, engineering_constants_reuss

def calculate_hybrid_properties(fiber_props: FiberProperties, matrix_props: MatrixProperties, fiber_volume_fraction: float) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    Vf = fiber_volume_fraction
    Vm = 1 - Vf

    ef1 = fiber_props.E1
    ef2 = fiber_props.E2
    ef3 = fiber_props.E2
    gf12 = fiber_props.G12
    gf13 = fiber_props.G12
    gf23 = fiber_props.E2 / (2 * (1 + fiber_props.nu23))
    nuf12 = fiber_props.nu12
    nuf13 = fiber_props.nu12
    nuf23 = fiber_props.nu23

    em1 = matrix_props.E1
    em2 = matrix_props.E1
    em3 = matrix_props.E1
    gm12 = matrix_props.E1 / (2 * (1 + matrix_props.nu))
    gm13 = gm12
    gm23 = gm12
    num12 = matrix_props.nu
    num13 = matrix_props.nu
    num23 = matrix_props.nu

    SHf_Temp = np.array([
        [ef1, nuf12, nuf13, 0, 0, 0],
        [-nuf12, 1 / ef2 - nuf12 * nuf12 / ef1, -nuf23 / ef2 - nuf13 * nuf13 / ef1, 0, 0, 0],
        [-nuf23, -nuf23 / ef2 - nuf12 * nuf12 / ef1, 1 / ef3 - nuf13 * nuf13 / ef1, 0, 0, 0],
        [0, 0, 0, 1 / gf23, 0, 0],
        [0, 0, 0, 0, 1 / gf13, 0],
        [0, 0, 0, 0, 0, 1 / gf12]
    ])

    SHm_Temp = np.array([
        [em1, num12, num13, 0, 0, 0],
        [-num12, 1 / em2 - num12 * num12 / em1, -num23 / em2 - num13 * num13 / em1, 0, 0, 0],
        [-num23, -num23 / em2 - num12 * num12 / em1, 1 / em3 - num13 * num13 / em1, 0, 0, 0],
        [0, 0, 0, 1 / gm23, 0, 0],
        [0, 0, 0, 0, 1 / gm13, 0],
        [0, 0, 0, 0, 0, 1 / gm12]
    ])

    SHs_Temp = SHf_Temp * Vf + SHm_Temp * Vm

    nu12 = SHs_Temp[0, 1]
    nu13 = SHs_Temp[0, 2]

    G12 = 1 / SHs_Temp[5, 5]
    G13 = 1 / SHs_Temp[4, 4]
    G23 = 1 / SHs_Temp[3, 3]

    E2 = 1 / (SHs_Temp[1, 1] + nu12 * nu12 / SHs_Temp[0, 0])
    E3 = 1 / (SHs_Temp[2, 2] + nu13 * nu13 / SHs_Temp[0, 0])
    nu23 = -E2 * (SHs_Temp[1, 2] + nu12 * nu12 / SHs_Temp[0, 0])

    hybrid_engineering_constants = {
        "E1": SHs_Temp[0, 0],
        "E2": E2,
        "E3": E3,
        "G12": G12,
        "G13": G13,
        "G23": G23,
        "nu12": nu12,
        "nu13": nu13,
        "nu23": nu23
    }

    Shs = np.array([
        [1 / SHs_Temp[0, 0], -nu12 / SHs_Temp[0, 0], -nu13 / SHs_Temp[0, 0], 0, 0, 0],
        [-nu12 / SHs_Temp[0, 0], 1 / E2, -nu23 / E2, 0, 0, 0],
        [-nu13 / SHs_Temp[0, 0], -nu23 / E2, 1 / E3, 0, 0, 0],
        [0, 0, 0, 1 / G23, 0, 0],
        [0, 0, 0, 0, 1 / G13, 0],
        [0, 0, 0, 0, 0, 1 / G12]
    ])

    Chs = np.linalg.inv(Shs)

    return Chs, Shs, hybrid_engineering_constants


# Request model for UDFRC properties calculation
class UDFRCPropertiesInput(BaseModel):
    fiber_E1: float = Field(..., gt=0, description="Longitudinal modulus of the fiber (must be greater than 0).")
    fiber_E2: float = Field(..., gt=0, description="Transverse modulus of the fiber (must be greater than 0).")
    fiber_G12: float = Field(..., gt=0, description="Shear modulus of the fiber (must be greater than 0).")
    fiber_nu12: float = Field(..., ge=-1, le=0.5, description="Poisson's ratio in the fiber plane (must be between -1 and 0.5).")
    fiber_nu23: float = Field(..., ge=-1, le=0.5, description="Poisson's ratio in the fiber transverse plane (must be between -1 and 0.5).")
    matrix_E1: float = Field(..., gt=0, description="Modulus of the matrix (must be greater than 0).")
    matrix_nu: float = Field(..., ge=-1, le=0.5, description="Poisson's ratio of the matrix (must be between -1 and 0.5).")
    fiber_volume_fraction: float = Field(..., ge=0, le=1, description="Volume fraction of fiber in the composite (must be between 0 and 1).")

# Response model for UDFRC properties calculation
class UDFRCPropertiesOutput(BaseModel):
    Voigt_Rules_of_Mixture: ThreeDimensionalPropertiesOutput
    Reuss_Rules_of_Mixture: ThreeDimensionalPropertiesOutput
    Hybrid_Rules_of_Mixture: ThreeDimensionalPropertiesOutput

@router.post("/udfrc-properties", response_model=UDFRCPropertiesOutput)
async def calcluate_udfrc_properties(data: UDFRCPropertiesInput):
    """
    API endpoint to calculate the UDFRC (Unidirectional Fiber Reinforced Composite) properties.

    Request Body:
    - fiber_E1 (float): Longitudinal modulus of the fiber.
    - fiber_E2 (float): Transverse modulus of the fiber.
    - fiber_G12 (float): Shear modulus of the fiber.
    - fiber_nu12 (float): Poisson's ratio in the fiber plane.
    - fiber_nu23 (float): Poisson's ratio in the fiber transverse plane.
    - matrix_E1 (float): Modulus of the matrix.
    - matrix_nu (float): Poisson's ratio of the matrix.
    - fiber_volume_fraction (float): Volume fraction of fiber in the composite (must be between 0 and 1).

    Returns a JSON object containing the Voigt, Reuss, and Hybrid rules of mixture with effective stiffness and compliance matrices, as well as engineering constants.
    
    ### Example `curl` Request:
    ```bash
    curl -X 'POST' \\
    'http://127.0.0.1:8000/api/v1/udfrc-properties' \\
    -H 'accept: application/json' \\
    -H 'Content-Type: application/json' \\
    -d '{
        "fiber_E1": 150000,
        "fiber_E2": 10000,
        "fiber_G12": 5000,
        "fiber_nu12": 0.3,
        "fiber_nu23": 0.25,
        "matrix_E1": 3500,
        "matrix_nu": 0.35,
        "fiber_volume_fraction": 0.6
        }'
    """

    fiber_E1 = data.fiber_E1
    fiber_E2 = data.fiber_E2
    fiber_G12 = data.fiber_G12
    fiber_nu12 = data.fiber_nu12
    fiber_nu23 = data.fiber_nu23
    matrix_E1 = data.matrix_E1
    matrix_nu = data.matrix_nu
    fiber_volume_fraction = data.fiber_volume_fraction

    try:
        fiber_props = FiberProperties(E1=fiber_E1, E2=fiber_E2, G12=fiber_G12, nu12=fiber_nu12, nu23=fiber_nu23)
        matrix_props = MatrixProperties(E1=matrix_E1, nu=matrix_nu)
        
        C_voigt, S_voigt, voigt_constants = calculate_voigt_properties(fiber_props, matrix_props, fiber_volume_fraction)
        C_reuss, S_reuss, reuss_constants = calculate_reuss_properties(fiber_props, matrix_props, fiber_volume_fraction)
        C_hybrid, S_hybrid, hybrid_constants = calculate_hybrid_properties(fiber_props, matrix_props, fiber_volume_fraction)
        
        result = {
            "Voigt_Rules_of_Mixture": {
                "Effective_3D_Stiffness_Matrix": C_voigt.tolist(),
                "Effective_3D_Compliance_Matrix": S_voigt.tolist(),
                "Engineering_Constants": voigt_constants
            },
            "Reuss_Rules_of_Mixture": {
                "Effective_3D_Stiffness_Matrix": C_reuss.tolist(),
                "Effective_3D_Compliance_Matrix": S_reuss.tolist(),
                "Engineering_Constants": reuss_constants
            },
            "Hybrid_Rules_of_Mixture": {
                "Effective_3D_Stiffness_Matrix": C_hybrid.tolist(),
                "Effective_3D_Compliance_Matrix": S_hybrid.tolist(),
                "Engineering_Constants": hybrid_constants
            }
        }
        
        return result
    except HTTPException as e:
        # Reraise HTTPException so it returns the correct status code
        raise e
    except Exception as e:
        print(f"Error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))
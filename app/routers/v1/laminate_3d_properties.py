from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import List, Tuple, Dict
import numpy as np
from math import radians, sin, cos, pi
from app.services.v1.layup_service import parse_layup_sequence
from app.models.v1.laminate_3d_properties_input import Laminate3DPropertiesInput
from app.models.v1.three_dimensional_properties_output import ThreeDimensionalPropertiesOutput

router = APIRouter()

def calculate_3d_properties(E1: float, E2: float, G12: float, nu12: float, nu23: float, layup_sequence: str, layer_thickness: float) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    layups = parse_layup_sequence(layup_sequence)
    n_ply = len(layups)
    thickness = layer_thickness

    bzi = [(-(n_ply + 1) * thickness) / 2 + i * thickness for i in range(1, n_ply + 1)]

    C = np.zeros((6, 6))

    for i in range(n_ply):
        layup = layups[i]
        e1 = E1
        e2 = E2
        g12 = G12
        nu12 = nu12
        nu23 = nu23
        e3 = e2
        g13 = g12
        g23 = e2 / (2 * (1 + nu23))
        nu13 = nu12
        angle_radian = layup * pi / 180
        s = sin(angle_radian)
        c = cos(angle_radian)
        
        Sp = np.array([
            [1 / e1, -nu12 / e1, -nu13 / e1, 0, 0, 0],
            [-nu12 / e1, 1 / e2, -nu23 / e2, 0, 0, 0],
            [-nu13 / e1, -nu23 / e2, 1 / e3, 0, 0, 0],
            [0, 0, 0, 1 / g23, 0, 0],
            [0, 0, 0, 0, 1 / g13, 0],
            [0, 0, 0, 0, 0, 1 / g12]
        ])

        Cp = np.linalg.inv(Sp)

        Rsigma = np.array([
            [c * c, s * s, 0, 0, 0, -2 * s * c],
            [s * s, c * c, 0, 0, 0, 2 * s * c],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, s, 0],
            [0, 0, 0, -s, c, 0],
            [s * c, -s * c, 0, 0, 0, c * c - s * s]
        ])
        
        C_single = Rsigma @ Cp @ Rsigma.T
        C += C_single

    C = C * (1 / n_ply)
    S = np.linalg.inv(C)

    engineering_constants = {
        "E1": 1 / S[0, 0],
        "E2": 1 / S[1, 1],
        "E3": 1 / S[2, 2],
        "G12": 1 / S[5, 5],
        "G13": 1 / S[4, 4],
        "G23": 1 / S[3, 3],
        "nu12": -S[0, 1] / S[0, 0],
        "nu13": -S[0, 2] / S[0, 0],
        "nu23": -S[1, 2] / S[1, 1]
    }

    return C, S, engineering_constants

@router.post("/laminate-3d-properties", response_model=ThreeDimensionalPropertiesOutput)
async def calculate_laminate_3d_properties(data: Laminate3DPropertiesInput):
    """
    API endpoint to calculate the 3D laminate properties of a composite material.

    Request Body:
    - E1 (float): Longitudinal modulus of the laminate.
    - E2 (float): Transverse modulus of the laminate.
    - G12 (float): In-plane shear modulus of the laminate.
    - nu12 (float): Poisson's ratio in the laminate plane.
    - nu23 (float): Poisson's ratio in the transverse plane.
    - layup_sequence (str): Layup sequence of the laminate (e.g., "[30/45]2s").
    - layer_thickness (float): Thickness of each laminate layer.

    Returns A JSON object containing the effective stiffness and compliance matrices, and engineering constants.
    
    ### Example `curl` Request:
    ```bash
    curl -X 'POST' \\
    'http://127.0.0.1:8000/api/v1/laminate-3d-properties' \\
    -H 'accept: application/json' \\
    -H 'Content-Type: application/json' \\
    -d '{
        "E1": 150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "nu23": 0.25,
        "layup_sequence": "[30/45]2s",
        "layer_thickness": 0.125
    }'
    
    """
    
    E1 = data.E1
    E2 = data.E2
    G12 = data.G12
    nu12 = data.nu12
    nu23 = data.nu23
    layup_sequence = data.layup_sequence
    layer_thickness = data.layer_thickness
    
    try:
        C, S, engineering_constants = calculate_3d_properties(
            E1, E2, G12, nu12, nu23, layup_sequence, layer_thickness
        )

        result = {
            "Effective_3D_Stiffness_Matrix": C.tolist(),
            "Effective_3D_Compliance_Matrix": S.tolist(),
            "Engineering_Constants": engineering_constants
        }
        
        return result
    except HTTPException as e:
        # Reraise HTTPException for validation errors
        raise e
    except Exception as e:
        print(f"Error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))
from fastapi import APIRouter, HTTPException
import numpy as np
from typing import List
from pydantic import BaseModel, Field
from math import radians, sin, cos

router = APIRouter()

# Define Pydantic model for input data
class LaminaEngineeringConstantsInput(BaseModel):
    E1: float = Field(..., gt=0, description="Young's Modulus in the longitudinal direction (must be greater than 0)")
    E2: float = Field(..., gt=0, description="Young's Modulus in the transverse direction (must be greater than 0)")
    G12: float = Field(..., gt=0, description="Shear Modulus in the plane of the lamina (must be greater than 0)")
    nu12: float = Field(..., ge=-1, le=0.5, description="Poisson's ratio in the plane of the lamina (between -1 and 0.5)")
    layup_angle: float = Field(..., description="The angle of the fiber layup in degrees")


# Define Pydantic model for output data
class LaminaEngineeringConstantsOutput(BaseModel):
    E_1: float
    E_2: float
    G_12: float
    nu_12: float
    eta_1_12: float
    eta_2_12: float
    Q: List[List[float]] = Field(..., description="Transformed stiffness matrix Q (2D array of floats)")
    S: List[List[float]] = Field(..., description="Transformed compliance matrix S (2D array of floats)")


@router.post("/lamina-engineering-constants", response_model=LaminaEngineeringConstantsOutput)
async def calculate_lamina_engineering_constants(data: LaminaEngineeringConstantsInput):
    """
    Calculate the lamina engineering constants based on provided material properties.
    
    - **E1**: Young's Modulus in the longitudinal direction.
    - **E2**: Young's Modulus in the transverse direction.
    - **G12**: Shear Modulus in the plane of the lamina.
    - **nu12**: Poisson's ratio in the plane of the lamina.
    - **layup_angle**: The angle of the fiber layup in degrees.
    
    This endpoint calculates:
    - Transformed stiffness matrix (Q).
    - Transformed compliance matrix (S).
    - Elastic moduli in the X and Y directions.
    - Shear modulus.
    - Interaction coefficients.
    
    Returns a JSON object with the calculated constants.

    ### Example `curl` Request:
    ```bash
    curl -X 'POST' \\
      'http://127.0.0.1:8000/api/v1/lamina-engineering-constants' \\
      -H 'accept: application/json' \\
      -H 'Content-Type: application/json' \\
      -d '{
        "E1": 150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "layup_angle": 45
      }'
    ```
    """
    E1 = data.E1
    E2 = data.E2
    G12 = data.G12
    nu12 = data.nu12
    layup_angle = data.layup_angle

    angle_radian = radians(layup_angle)
    
    S = np.array([
        [1 / E1, -nu12 / E1, 0],
        [-nu12 / E1, 1 / E2, 0],
        [0, 0, 1 / G12]
    ])
    
    print("Matrix S:", S)
    
    try:
        C = np.linalg.inv(S)
    except np.linalg.LinAlgError as e:
        print("Matrix S is singular and cannot be inverted:", e)
        raise HTTPException(status_code=500, detail="Matrix S is singular and cannot be inverted")
    
    s = sin(angle_radian)
    c = cos(angle_radian)
    
    T_epsilon = np.array([
        [c * c, s * s, -s * c],
        [s * s, c * c, s * c],
        [2 * s * c, -2 * s * c, c * c - s * s]
    ])
    
    T_sigma = np.array([
        [c * c, s * s, -2 * s * c],
        [s * s, c * c, 2 * s * c],
        [s * c, -s * c, c * c - s * s]
    ])
    
    C_bar = np.dot(np.dot(T_sigma.T, C), T_sigma)
    S_bar = np.dot(np.dot(T_epsilon.T, S), T_epsilon)
    
    E_x = 1 / S_bar[0, 0]
    E_y = 1 / S_bar[1, 1]
    G_xy = 1 / S_bar[2, 2]
    nu_xy = -S_bar[0, 1] * E_x
    eta_x_xy = S_bar[2, 0] * E_x
    eta_y_xy = S_bar[2, 1] * E_y
    
    result = {
        "E_1": E_x,
        "E_2": E_y,
        "G_12": G_xy,
        "nu_12": nu_xy,
        "eta_1_12": eta_x_xy,
        "eta_2_12": eta_y_xy,
        "Q": C_bar.tolist(),
        "S": S_bar.tolist()
    }
    
    return result
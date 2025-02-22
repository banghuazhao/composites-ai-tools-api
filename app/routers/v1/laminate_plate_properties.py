from fastapi import APIRouter, HTTPException
from typing import List, Tuple
import numpy as np
from math import radians, sin, cos, pi
from pydantic import BaseModel, Field
from app.services.v1.layup_service import parse_layup_sequence
from app.models.v1.laminate_plate_properties_input import LaminatePlatePropertiesInput

router = APIRouter()

def calculate_laminate_properties(E1: float, E2: float, G12: float, nu12: float, layup_sequence: str, layer_thickness: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    layups = parse_layup_sequence(layup_sequence)
    n_ply = len(layups)
    thickness = layer_thickness

    bzi = [(-(n_ply + 1) * thickness) / 2 + i * thickness for i in range(1, n_ply + 1)]

    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    for i in range(n_ply):
        layup = layups[i]
        angle_radian = layup * pi / 180
        s = sin(angle_radian)
        c = cos(angle_radian)

        Sep = np.array([
            [1 / E1, -nu12 / E1, 0],
            [-nu12 / E1, 1 / E2, 0],
            [0, 0, 1 / G12]
        ])

        Qep = np.linalg.inv(Sep)

        Rsigmae = np.array([
            [c * c, s * s, -2 * s * c],
            [s * s, c * c, 2 * s * c],
            [s * c, -s * c, c * c - s * s]
        ])

        Qe = np.dot(np.dot(Rsigmae, Qep), Rsigmae.T)

        A += Qe * thickness
        B += Qe * thickness * bzi[i]
        D += Qe * (thickness * bzi[i] * bzi[i] + (thickness ** 3) / 12)

    h = n_ply * thickness

    try:
        Ses = np.linalg.inv(A) * h
        Sesf = np.linalg.inv(D) * (h ** 3 / 12)
    except np.linalg.LinAlgError:
        raise HTTPException(status_code=500, detail="Singular matrix encountered during inversion")

    in_plane_properties = {
        "e1": 1 / Ses[0, 0],
        "e2": 1 / Ses[1, 1],
        "g12": 1 / Ses[2, 2],
        "nu12": -Ses[0, 1] / Ses[0, 0],
        "eta121": -Ses[0, 2] / Ses[2, 2],
        "eta122": -Ses[1, 2] / Ses[2, 2]
    }

    flexural_properties = {
        "e1": 1 / Sesf[0, 0],
        "e2": 1 / Sesf[1, 1],
        "g12": 1 / Sesf[2, 2],
        "nu12": -Sesf[0, 1] / Sesf[0, 0],
        "eta121": -Sesf[0, 2] / Sesf[2, 2],
        "eta122": -Sesf[1, 2] / Sesf[2, 2]
    }

    return A, B, D, in_plane_properties, flexural_properties

class LaminateProperties(BaseModel):
    e1: float = Field(..., description="In-plane or flexural modulus in the longitudinal direction.")
    e2: float = Field(..., description="In-plane or flexural modulus in the transverse direction.")
    g12: float = Field(..., description="Shear modulus in the 1-2 plane.")
    nu12: float = Field(..., description="Poisson's ratio in the 1-2 plane.")
    eta121: float = Field(..., description="Shear coupling term (eta121).")
    eta122: float = Field(..., description="Shear coupling term (eta122).")

class LaminatePlatePropertiesResponse(BaseModel):
    A: list = Field(..., description="In-plane stiffness matrix (6x6).")
    B: list = Field(..., description="Coupling stiffness matrix (6x6).")
    D: list = Field(..., description="Flexural stiffness matrix (6x6).")
    in_plane_properties: LaminateProperties = Field(..., description="In-plane engineering properties.")
    flexural_properties: LaminateProperties = Field(..., description="Flexural engineering properties.")


@router.post("/laminate-plate-properties", response_model=LaminatePlatePropertiesResponse)
async def calculate_laminate_plate_properties(data: LaminatePlatePropertiesInput):
    """
    API endpoint to calculate the laminate plate properties of a composite material.
    
    This endpoint accepts material properties and a layup sequence as input, and returns the laminate stiffness matrices
    and other related properties such as in-plane and flexural properties.

    Payload: The input model containing:
    E1 (float): Young's modulus in the longitudinal direction.
    E2 (float): Young's modulus in the transverse direction.
    G12 (float): Shear modulus in the plane of the laminate.
    nu12 (float): Poisson's ratio in the plane of the laminate.
    layup_sequence (str): Layup sequence of the laminate, e.g., "[45/90/-45]s".
    layer_thickness (float): Thickness of each layer in the laminate.
    
    Returns a dictionary containing:
    - "A" (List[List[float]]): The in-plane stiffness matrix (6x6).
    - "B" (List[List[float]]): The coupling stiffness matrix (6x6).
    - "D" (List[List[float]]): The flexural stiffness matrix (6x6).
    - "in_plane_properties" (dict): Dictionary containing the in-plane engineering properties.
    - "flexural_properties" (dict): Dictionary containing the flexural engineering properties.
    
    ### Example `curl` Request:
    ```bash
    curl -X 'POST' \\
    'http://127.0.0.1:8000/api/v1/laminate-plate-properties' \\
    -H 'accept: application/json' \\
    -H 'Content-Type: application/json' \\
    -d '{
        "E1": 150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "layup_sequence": "[45/90/-45]s",
        "layer_thickness": 0.125
        }'
    """

    E1 = data.E1
    E2 = data.E2
    G12 = data.G12
    nu12 = data.nu12
    layup_sequence = data.layup_sequence
    layer_thickness = data.layer_thickness
    
    try:
        A, B, D, in_plane_properties, flexural_properties = calculate_laminate_properties(E1, E2, G12, nu12, layup_sequence, layer_thickness)
        
        result = {
            "A": A.tolist(),
            "B": B.tolist(),
            "D": D.tolist(),
            "in_plane_properties": in_plane_properties,
            "flexural_properties": flexural_properties
        }

        return result
    except HTTPException as e:
        # Reraise HTTPException so it returns the correct status code
        raise e
    except Exception as e:
        print(f"Error: {str(e)}")
        raise HTTPException(status_code=500, detail="An internal server error occurred: " + str(e))
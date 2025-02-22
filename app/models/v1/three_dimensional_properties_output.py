from pydantic import BaseModel, Field
from typing import List

# Model for the engineering constants
class EngineeringConstants(BaseModel):
    E1: float
    E2: float
    E3: float
    G12: float
    G13: float
    G23: float
    nu12: float
    nu13: float
    nu23: float

# Model for Voigt, Reuss, Hybrid rules of mixture
class ThreeDimensionalPropertiesOutput(BaseModel):
    Effective_3D_Stiffness_Matrix: List[List[float]] = Field(..., description="Effective 3D stiffness matrix (6x6).")
    Effective_3D_Compliance_Matrix: List[List[float]] = Field(..., description="Effective 3D compliance matrix (6x6).")
    Engineering_Constants: EngineeringConstants = Field(..., description="Engineering constants.")

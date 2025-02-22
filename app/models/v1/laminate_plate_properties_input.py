from pydantic import BaseModel, Field

# Request Model
class LaminatePlatePropertiesInput(BaseModel):
    E1: float = Field(..., gt=0, description="Longitudinal modulus of the laminate (must be greater than 0).")
    E2: float = Field(..., gt=0, description="Transverse modulus of the laminate (must be greater than 0).")
    G12: float = Field(..., gt=0, description="In-plane shear modulus of the laminate (must be greater than 0).")
    nu12: float = Field(..., ge=-1, le=0.5, description="Poisson's ratio in the laminate plane (must be between -1 and 0.5).")
    layup_sequence: str = Field(..., description="Layup sequence of the laminate (e.g., '[30/45]2s').")
    layer_thickness: float = Field(..., gt=0, description="Thickness of each laminate layer (must be greater than 0).")

from fastapi.testclient import TestClient
from app.main import app  # Adjust this import based on your project structure

client = TestClient(app)

def test_calculate_laminate_3d_properties_success():
    # Valid input data
    payload = {
        "E1": 150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "nu23": 0.25,
        "layup_sequence": "[0/90/45]",
        "layer_thickness": 0.125
    }

    # Send a POST request
    response = client.post("/api/v1/laminate-3d-properties", json=payload)

    # Assert the status code is 200 (OK)
    assert response.status_code == 200

    # Check the response content
    response_data = response.json()
    assert "Effective_3D_Stiffness_Matrix" in response_data
    assert "Effective_3D_Compliance_Matrix" in response_data
    assert "Engineering_Constants" in response_data


def test_calculate_laminate_3d_properties_invalid_data():
    # Test with invalid input data (negative E1 value)
    payload = {
        "E1": -150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "nu23": 0.25,
        "layup_sequence": "0/90/45",
        "layer_thickness": 0.125
    }

    # Send a POST request
    response = client.post("/api/v1/laminate-3d-properties", json=payload)

    assert response.status_code == 422

def test_calculate_laminate_3d_properties_missing_param():
    # Test with missing parameter (no 'G12')
    payload = {
        "E1": 150000,
        "E2": 10000,
        "nu12": 0.3,
        "nu23": 0.25,
        "layup_sequence": "0/90/45",
        "layer_thickness": 0.125
    }

    # Send a POST request
    response = client.post("/api/v1/laminate-3d-properties", json=payload)

    assert response.status_code == 422

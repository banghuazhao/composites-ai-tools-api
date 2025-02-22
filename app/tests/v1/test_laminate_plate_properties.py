from fastapi.testclient import TestClient
from app.main import app  # Replace with the actual name of your main FastAPI app file

# Create a TestClient using the FastAPI app instance
client = TestClient(app)

def test_get_laminate_properties_success():
    # Sample input data
    payload = {
        "E1": 150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "layup_sequence": "0/90/45",
        "layer_thickness": 0.125
    }

    # Send a POST request
    response = client.post("/api/v1/laminate-plate-properties", json=payload)

    # Assert the status code is 200 (OK)
    assert response.status_code == 200

    # Check the response JSON structure
    response_data = response.json()
    assert "A" in response_data
    assert "B" in response_data
    assert "D" in response_data
    assert "in_plane_properties" in response_data
    assert "flexural_properties" in response_data

def test_get_laminate_properties_invalid_data():
    # Test with invalid input data (negative E1 value)
    payload = {
        "E1": -150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "layup_sequence": "0/90/45",
        "layer_thickness": 0.125
    }

    # Send a POST request
    response = client.post("/api/v1/laminate-plate-properties", json=payload)

    assert response.status_code == 422

def test_get_laminate_properties_missing_param():
    # Test with missing parameter (no 'G12')
    payload = {
        "E1": 150000,
        "E2": 10000,
        "nu12": 0.3,
        "layup_sequence": "0/90/45",
        "layer_thickness": 0.125
    }

    # Send a POST request
    response = client.post("/api/v1/laminate-plate-properties", json=payload)

    assert response.status_code == 422
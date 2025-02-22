from fastapi.testclient import TestClient
from app.main import app  # Adjust this import according to your project structure

client = TestClient(app)

def test_calculate_lamina_constants_success():
    # Valid input data
    payload = {
        "E1": 150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "layup_angle": 45
    }

    # Send a POST request
    response = client.post("/api/v1/lamina-engineering-constants", json=payload)

    # Assert the status code is 200 (OK)
    assert response.status_code == 200

    # Check the response content
    response_data = response.json()
    assert "E_1" in response_data
    assert "E_2" in response_data
    assert "G_12" in response_data
    assert "nu_12" in response_data
    assert "eta_1_12" in response_data
    assert "eta_2_12" in response_data
    assert "Q" in response_data
    assert "S" in response_data


def test_calculate_lamina_constants_invalid_data():
    # Test with invalid input data (negative E1 value)
    payload = {
        "E1": -150000,
        "E2": 10000,
        "G12": 5000,
        "nu12": 0.3,
        "layup_angle": 45
    }

    # Send a POST request
    response = client.post("/api/v1/lamina-engineering-constants", json=payload)

    # Assert the status code is 422 (Unprocessable Entity)
    assert response.status_code == 422

def test_calculate_lamina_constants_missing_param():
    # Test with missing parameter (no 'G12')
    payload = {
        "E1": 150000,
        "E2": 10000,
        "nu12": 0.3,
        "layup_angle": 45
    }

    # Send a POST request
    response = client.post("/api/v1/lamina-engineering-constants", json=payload)

    # Assert the status code is 422 (Unprocessable Entity)
    assert response.status_code == 422

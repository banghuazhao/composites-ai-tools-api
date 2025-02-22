from fastapi.testclient import TestClient
from app.main import app  # Adjust this import based on your project structure

client = TestClient(app)

def test_udfrc_properties_success():
    # Valid input data
    payload = {
        "fiber_E1": 150000,
        "fiber_E2": 10000,
        "fiber_G12": 5000,
        "fiber_nu12": 0.3,
        "fiber_nu23": 0.25,
        "matrix_E1": 2500,
        "matrix_nu": 0.35,
        "fiber_volume_fraction": 0.6
    }

    # Send a POST request
    response = client.post("/api/v1/udfrc-properties", json=payload)

    # Assert the status code is 200 (OK)
    assert response.status_code == 200

    # Check the response content
    response_data = response.json()
    assert "Voigt_Rules_of_Mixture" in response_data
    assert "Reuss_Rules_of_Mixture" in response_data
    assert "Hybrid_Rules_of_Mixture" in response_data

def test_udfrc_properties_invalid_modulus():
    # Invalid input: negative modulus
    payload = {
        "fiber_E1": -150000,  # Negative value
        "fiber_E2": 10000,
        "fiber_G12": 5000,
        "fiber_nu12": 0.3,
        "fiber_nu23": 0.25,
        "matrix_E1": 2500,
        "matrix_nu": 0.35,
        "fiber_volume_fraction": 0.6
    }

    # Send a POST request
    response = client.post("/api/v1/udfrc-properties", json=payload)

    assert response.status_code == 422

def test_udfrc_properties_invalid_poisson_ratio():
    # Invalid input: Poisson's ratio out of range
    payload = {
        "fiber_E1": 150000,
        "fiber_E2": 10000,
        "fiber_G12": 5000,
        "fiber_nu12": 1.0,  # Invalid Poisson's ratio
        "fiber_nu23": 0.25,
        "matrix_E1": 2500,
        "matrix_nu": 0.35,
        "fiber_volume_fraction": 0.6
    }

    # Send a POST request
    response = client.post("/api/v1/udfrc-properties", json=payload)

    assert response.status_code == 422

def test_udfrc_properties_invalid_fiber_volume_fraction():
    # Invalid input: fiber volume fraction out of range
    payload = {
        "fiber_E1": 150000,
        "fiber_E2": 10000,
        "fiber_G12": 5000,
        "fiber_nu12": 0.3,
        "fiber_nu23": 0.25,
        "matrix_E1": 2500,
        "matrix_nu": 0.35,
        "fiber_volume_fraction": 1.5  # Invalid fraction (> 1)
    }

    # Send a POST request
    response = client.post("/api/v1/udfrc-properties", json=payload)

    assert response.status_code == 422

def test_udfrc_properties_missing_parameters():
    # Missing 'fiber_G12' parameter
    payload = {
        "fiber_E1": 150000,
        "fiber_E2": 10000,
        "fiber_nu12": 0.3,
        "fiber_nu23": 0.25,
        "matrix_E1": 2500,
        "matrix_nu": 0.35,
        "fiber_volume_fraction": 0.6
    }

    # Send a POST request
    response = client.post("/api/v1/udfrc-properties", json=payload)

    assert response.status_code == 422
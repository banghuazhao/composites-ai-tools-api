# Composites-AI-Tools-API
[![Version](https://img.shields.io/github/v/release/banghuazhao/composites-ai-tools-api)](https://github.com/banghuazhao/composites-ai-tools-api/releases)
[![Build Status](https://github.com/banghuazhao/composites-ai-tools-api/actions/workflows/ci.yml/badge.svg)](https://github.com/banghuazhao/composites-ai-tools-api/actions)
[![codecov](https://codecov.io/gh/banghuazhao/composites-ai-tools-api/branch/main/graph/badge.svg)](https://codecov.io/gh/banghuazhao/composites-ai-tools-api)
[![Last Commit](https://img.shields.io/github/last-commit/banghuazhao/composites-ai-tools-api)](https://github.com/banghuazhao/composites-ai-tools-api/composites-ai-tools/main)


**Composites-AI-Tools-API** is a FastAPI-based backend designed to perform calculations related to composite laminate materials. It exposes various endpoints for computing engineering constants, plate properties, 3D laminate properties, and more. This API is suitable for researchers, engineers, and developers working in composite materials and structural analysis.

## Table of Contents

- [Features](#features)
- [Technologies](#technologies)
- [Installation Locally](#installation-locally)
- [API Documentation](#api-documentation)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

---

## Features

- **Versioned API:** Support for multiple versions of the API, with versioning (v1, v2) for backward compatibility.
- **Material Property Calculations:** Endpoints to calculate various laminate engineering constants and composite properties.
- **Lightweight & Fast:** Built using FastAPI, ensuring high performance and minimal latency.
- **Modular Architecture:** Separation of concerns using routers, services, and utility functions for better maintainability.
- **Extensible:** Easy to add new endpoints and services for additional material computations.

---

## Technologies

- **FastAPI**: A modern, fast (high-performance) web framework for building APIs with Python 3.10.12+ based on standard Python-type hints.
- **Python 3.10.12**: Programming language used to implement the logic.
- **Pytest**: A robust testing framework for Python.
- **Uvicorn**: ASGI server used for serving FastAPI applications.
- **Virtual Environment (venv)**: For dependency isolation.

---

## Installation Locally

To set up the project locally, follow these steps:

### Prerequisites
- Python 3.10.12 or higher.
- Git installed on your system.

### Step-by-Step Guide

1. **Clone the repository**:
```bash
git clone https://github.com/banghuazhao/composites-ai-tools-api.git
cd composites-ai-tools-api
```

2. **Create and activate a virtual environment**:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install the dependencies**:
```bash
pip install -r requirements.txt
```

4. **Run the FastAPI application locally**:
```bash
uvicorn app.main:app --reload
```

This will start the Composites AI Tools API server at http://127.0.0.1:8000.

---

## Redis Setup

APIs may use Redis for data storage. In the code of an API that uses Redis, the Redis connection can be configured via an environment variable:

```python
redis_url = os.getenv("REDIS_URL")
```

To run the backend with your own Redis server, set the REDIS_URL environment variable:

```bash
export REDIS_URL=redis://localhost:6379
```

For production deployments, define this variable in your `systemd` service file, `.env` file, or deployment configuration so the environment variable `REDIS_URL` loads automatically when the server starts.

## API Documentation
FastAPI automatically generates interactive API documentation using **Swagger** and **Redoc**. Once the server is running, you can access the documentation at:

* Swagger UI: http://127.0.0.1:8000/docs
* ReDoc: http://127.0.0.1:8000/redoc

## Sample Endpoints
* **GET /**: Root endpoint, returns a welcome message.
* **POST /api/v1/lamina_engineering_constants**: Calculates the lamina engineering constants based on the provided material properties.
* **POST /api/v1/laminate_plate_properties**: Calculates the laminate plate properties.
* **POST /api/v1/laminate_3d_properties**: Computes the 3D properties of laminate.
* **POST /api/v1/udfrc_properties**: Calculates UDFRC properties.

## Testing
The repository includes unit tests using `pytest`. To run the tests:

1. Ensure you are in the virtual environment.
2. Run the following command:
```bash
PYTHONPATH=./ pytest
```

## Contributing
We welcome contributions to **Composites-AI-Tools-API**! Here's how you can help:

1. **Fork the repository**.
2. **Create a new feature branch**:
```bash
git checkout -b feature/your-feature-name
```
3. **Create your function**
Create a single file in in `app/routers/v1` folder.
You can refer to `lamina_engineering_constants.py` as an example
Please define input and output model and define your function like this:
```python
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
```
3. **Commit your changes**:
```bash
git commit -m "Add your message"
```
4. **Push to the branch**:
```bash
git push origin feature/your-feature-name
```
5. **Create a pull request**: Describe the changes you’ve made.


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

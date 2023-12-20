"""Settings for PK models."""
from pathlib import Path

BASE_DIR = Path(__file__).parent
MODELS_DIR = BASE_DIR / "models" / "results"

DATA_DIR = BASE_DIR.parent / "measurement_files"
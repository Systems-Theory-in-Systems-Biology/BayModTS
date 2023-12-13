"""Module for creating all PK models"""

from pk_models import MODELS_DIR
from pk_models.models import simple_pk
from pk_models.models.templates import create_pk_model

if __name__ == "__main__":


    # FIXME: generalize for all models

    # test model
    create_pk_model(model=simple_pk._m, models_dir=MODELS_DIR, visualize=True)

    # caffeine
    # TODO

    # midazolam
    # TODO

    # codeine
    # TODO

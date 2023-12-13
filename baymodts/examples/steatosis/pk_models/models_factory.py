"""Module for creating all PK models"""
from pk_models import MODELS_DIR
from pk_models.models import simple_pk

if __name__ == "__main__":


    # FIXME: generalize for all models

    # test model
    simple_pk.create_simple_pk(models_dir=MODELS_DIR, visualize=True)

    # caffeine
    # TODO

    # midazolam
    # TODO

    # codeine
    # TODO

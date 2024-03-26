"""AMICI-generated module for model Control_rtf"""

from pathlib import Path

import amici

# Ensure we are binary-compatible, see #556
if "0.21.2" != amici.__version__:
    raise amici.AmiciVersionError(
        f"Cannot use model `Control_rtf` in {Path(__file__).parent}, "
        "generated with amici==0.21.2, "
        f"together with amici=={amici.__version__} "
        "which is currently installed. To use this model, install "
        "amici==0.21.2 or re-import the model with the amici "
        "version currently installed."
    )

from .Control_rtf import *  # noqa: F403, F401
from .Control_rtf import getModel as get_model  # noqa: F401

__version__ = "0.1.0"

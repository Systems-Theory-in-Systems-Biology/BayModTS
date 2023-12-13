"""Example simulations for caffeine model."""
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import roadrunner
from matplotlib import pyplot as plt
from sbmlutils.console import console

labels = {
    "caf_gut": "caffeine (gut)",
    "caf_cent": "caffeine (central)",
    "caf_peri": "caffeine (peripheral)",
    "caf_plasma": "caffeine (plasma)",
}


def reference_simulation(r: roadrunner.RoadRunner) -> None:
    """Reference simulation."""

    console.rule("Caffeine simulation", style="white")
    panel_size = 3
    f, (ax1, ax2) = plt.subplots(
        nrows=1, ncols=2, figsize=(2*panel_size, panel_size), dpi=300,
        layout="constrained"
    )

    # simulation
    r.resetAll()
    tend
    s = r.simulate(start=0, end=10, steps=400)  # [min]
    df: pd.DataFrame = pd.DataFrame(s, columns=s.colnames)

    for sid in ["caf_gut", "caf_cent", "caf_peri"]:
        ax1.plot(
            df.time/60,  # [min] -> [hr]
            df[f"[{sid}]"],
            label=labels[sid],
        )

        # ax.legend()
        ax1.set_title("Reference simulation", fontweight="bold", fontsize="15")
        ax1.set_xlabel("time [hr]", fontweight="bold", fontsize="12")
        ax1.set_ylabel("concentration [mmole/l]", fontweight="bold", fontsize="12")

    sid = "caf_plasma"
    ax2.plot(
        df.time / 60,  # [min] -> [hr]
        df[sid],
        label=labels[sid],
    )
    ax2.set_title("Reference observables", fontweight="bold", fontsize="15")
    ax2.set_xlabel("time [hr]", fontweight="bold", fontsize="12")
    ax2.set_ylabel("caffeine [ng/ml]", fontweight="bold", fontsize="12")

    for ax in (ax1, ax2):
        ax.legend()
    plt.show()


if __name__ == "__main__":
    # load_model
    from pk_models import MODELS_DIR
    model_path: Path = MODELS_DIR / "caffeine_pk.xml"
    r: roadrunner.RoadRunner = roadrunner.RoadRunner(str(model_path))
    model: roadrunner.ExecutableModel = r.model

    # selections
    selections = ["time"] + [f"[{sid}]" for sid in model.getFloatingSpeciesIds()] + [sid for sid in model.getGlobalParameterIds()]
    console.print(selections)
    r.selections = selections

    # changes
    changes: Dict[str, float] = {
        "k": 1E-2,  # [l/min],
        "Q": 1,  # [l/min],
        "CL": 1E-2,  # [l/min],
    }

    reference_simulation(r)

"""Example simulations for caffeine model."""
from pathlib import Path

import numpy as np
import pandas as pd
import roadrunner
from matplotlib import pyplot as plt
from sbmlutils.console import console

labels = {
    "caf_gut": "caffeine (gut)",
    "caf_cent": "caffeine (plasma)",
    "caf_peri": "caffeine (peripheral)",
}

def reference_simulation(r: roadrunner.RoadRunner) -> None:
    """Reference simulation."""

    console.rule("Caffeine simulation", style="white")
    f, ax = plt.subplots(
        nrows=1, ncols=1, figsize=(5, 5), dpi=300,
        layout="constrained"
    )

    # simulation
    r.resetAll()
    s = r.simulate(start=0, end=10, steps=400)  # [min]
    df: pd.DataFrame = pd.DataFrame(s, columns=s.colnames)

    for sid in ["caf_gut", "caf_cent", "caf_peri"]:
        linewidth = 1.0
        if sid == "caf_cent":
            linewidth = 2.0

        ax.plot(
            df.time/60,  # [min] -> [hr]
            df[f"[{sid}]"],
            label=labels[sid],
            linewidth=linewidth,
        )

        # ax.legend()
        ax.set_title("Reference simulation", fontweight="bold", fontsize="15")
        ax.set_xlabel("time [min]", fontweight="bold", fontsize="12")
        ax.set_ylabel("concentration [mmole/l]", fontweight="bold", fontsize="12")

    ax.legend()
    plt.show()


if __name__ == "__main__":
    # load_model
    from pk_models import MODELS_DIR
    model_path: Path = MODELS_DIR / "caffeine_pk.xml"
    r: roadrunner.RoadRunner = roadrunner.RoadRunner(str(model_path))
    reference_simulation(r)

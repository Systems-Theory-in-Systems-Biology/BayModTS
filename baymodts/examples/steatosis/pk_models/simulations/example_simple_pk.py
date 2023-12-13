"""Example simulations and parameter scans with simple PK model."""
from pathlib import Path

import numpy as np
import pandas as pd
import roadrunner
from matplotlib import pyplot as plt
from sbmlutils.console import console


def reference_simulation(r: roadrunner.RoadRunner) -> None:
    """Reference simulation."""

    console.rule("Simple PK simulation", style="white")
    f, ax = plt.subplots(
        nrows=1, ncols=1, figsize=(5, 5), dpi=300,
        layout="constrained"
    )

    # simulation
    r.resetAll()
    s = r.simulate(start=0, end=10, steps=400)  # [min]
    df: pd.DataFrame = pd.DataFrame(s, columns=s.colnames)

    for sid in ["[y_gut]", "[y_cent]", "[y_peri]"]:
        linewidth = 1.0
        if sid == "[y_cent]":
            linewidth = 2.0

        ax.plot(
            df.time, df[sid],
            label=sid,
            linewidth=linewidth,
        )

        # ax.legend()
        ax.set_title("Reference simulation", fontweight="bold", fontsize="15")
        ax.set_xlabel("time [min]", fontweight="bold", fontsize="12")
        ax.set_ylabel("concentration [mM]", fontweight="bold", fontsize="12")

    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # load_model
    from pk_models import MODELS_DIR
    simple_pk_path: Path = MODELS_DIR / "simple_pk.xml"
    r: roadrunner.RoadRunner = roadrunner.RoadRunner(str(simple_pk_path))
    reference_simulation(r)

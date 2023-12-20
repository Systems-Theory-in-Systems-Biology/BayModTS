"""Example simulations for PK model."""
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import roadrunner
from matplotlib import pyplot as plt
from sbmlutils.console import console

from pk_models import DATA_DIR

labels = {
    "caf_gut": "caffeine (gut)",
    "caf_cent": "caffeine (central)",
    "caf_peri": "caffeine (peripheral)",
    "caf_plasma": "caffeine (plasma)",

    "mid_gut": "midazolam (gut)",
    "mid_cent": "midazolam (central)",
    "mid_peri": "midazolam (peripheral)",
    "mid_plasma": "midazolam (plasma)",

    "cod_gut": "codeine (gut)",
    "cod_cent": "codeine (central)",
    "cod_peri": "codeine (peripheral)",
    "cod_plasma": "codeine (plasma)",
}

conditions = ["Control", "2Wks", "4Wks"]
condition_colors = {
    "Control": "black",
    "2Wks": "tab:orange",
    "4Wks": "tab:red",
}


def reference_simulation(sid: str, name: str) -> None:
    """Reference simulation."""
    console.rule(f"{name} simulation", style="white")

    # load_model
    from pk_models import MODELS_DIR
    model_path: Path = MODELS_DIR / f"{name}_pk.xml"
    r: roadrunner.RoadRunner = roadrunner.RoadRunner(str(model_path))
    model: roadrunner.ExecutableModel = r.model

    # selections
    selections = ["time"] + [f"[{sid}]" for sid in model.getFloatingSpeciesIds()] + [sid for sid in model.getGlobalParameterIds()]
    console.print(selections)
    r.selections = selections

    # changes
    changes: Dict[str, float] = {
        # "kabs": 1E-2,  # [l/hr],
        # "Q": 1E-1,  # [l/hr],
        # "CL": 1E-2,  # [l/hr],
    }

    panel_size = 4
    f, (ax1, ax2) = plt.subplots(
        nrows=1, ncols=2, figsize=(2*panel_size, panel_size), dpi=300,
        layout="constrained"
    )

    # simulation
    r.resetAll()
    for key, value in changes.items():
        r.setValue(key, value)

    s = r.simulate(start=0, end=7, steps=400)  # [min]
    df: pd.DataFrame = pd.DataFrame(s, columns=s.colnames)

    for key in [f"{sid}_gut", f"{sid}_cent", f"{sid}_peri"]:
        ax1.plot(
            df.time,
            df[f"[{key}]"],
            label=labels[key],
        )

        # ax.legend()
        ax1.set_title("Reference simulation", fontweight="bold", fontsize="15")
        ax1.set_xlabel("time [hr]", fontweight="bold", fontsize="12")
        ax1.set_ylabel("concentration [mmole/l]", fontweight="bold", fontsize="12")

    key = f"{sid}_plasma"
    ax2.plot(
        df.time,
        df[key],
        label=labels[key],
        color="black",
    )
    ax2.set_title("Reference observables", fontweight="bold", fontsize="15")
    ax2.set_xlabel("time [hr]", fontweight="bold", fontsize="12")
    ax2.set_ylabel(f"{name} [ng/ml]", fontweight="bold", fontsize="12")

    # load data and plot data
    conditions = ["Control", "2Wks", "4Wks"]
    substance = f"{name.title()}"
    for condition in conditions:
        csv_path = DATA_DIR / f"{substance}_{condition}_measurement_table.tsv"
        df = pd.read_csv(csv_path, sep="\t")
        # console.print(df)

        # mean +- SD
        times = []
        means = []
        stds = []

        for time, df_time in df.groupby("time"):
            times.append(df_time.time.mean())
            means.append(df_time.measurement.mean())
            stds.append(df_time.measurement.std())

        ax2.errorbar(
            x=times,
            y=means,
            yerr=stds,
            marker="s",
            color=condition_colors[condition],
            markeredgecolor="black",
            linestyle="-",
            label=condition,
            # markersize=5,
        )

    for condition in conditions:
        # individual data points
        ax2.plot(
            df.time, df.measurement,
            marker="o",
            color=condition_colors[condition],
            markeredgecolor="black",
            linestyle="",
            label=condition,
            markersize=2,
            alpha=0.8
        )

    for ax in (ax1, ax2):
        ax.legend()

    plt.show()
    f.savefig(MODELS_DIR / f"{name}_example.png", bbox_inches="tight")


if __name__ == "__main__":
    reference_simulation(sid="caf", name="caffeine")
    reference_simulation(sid="mid", name="midazolam")
    reference_simulation(sid="cod", name="codeine")

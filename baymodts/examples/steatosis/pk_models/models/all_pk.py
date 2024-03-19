"""Caffeine, codeine, midazolam PK model."""
from dataclasses import dataclass
from pathlib import Path
import os

from sbmlutils.console import console
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.examples.templates import terms_of_use
from sbmlutils.factory import *
from sbmlutils.metadata import *
import annotations

from templates import create_pk_model


class U(Units):
    """UnitDefinitions."""

    mmole = UnitDefinition("mmole")
    hr = UnitDefinition("hr")
    mg = UnitDefinition("mg")
    mug = UnitDefinition("mug", "µg")
    per_hr = UnitDefinition("per_hr", "1/hr")
    per_hr_l = UnitDefinition("per_hr_l", "1/hr/liter")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_hr = UnitDefinition("mmole_per_hr", "mmole/hr")
    mmole_per_hr_l = UnitDefinition("mmole_per_hr_l", "mmole/hr/liter")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")
    l_per_hr = UnitDefinition("l_per_hr", "l/hr")
    l_per_hr_mmole = UnitDefinition("l_per_hr_mmole", "l/hr/mmole")
    mM = UnitDefinition("mM", "mmole/liter")
    ng_per_ml = UnitDefinition("ng_per_ml", "mg/ml")
    ng_per_mug = UnitDefinition("ng_per_mug", "ng/µg")


def create_model_from_info(
        sid: str, name: str, Mr: float, dose_bw: float
):
    bodyweight = 0.029  # kg

    _m = Model(
        sid=f"{name}_pk",
        name=f"{name.title()} pharmacokinetics model",
        notes=f"""
        # Model of absorption and distribution of {name}.
        
        Simple two-compartment pharmacokinetics model for {name}.
        """
        + terms_of_use,
        creators=[
            Creator(
                familyName="König",
                givenName="Matthias",
                email="koenigmx@hu-berlin.de",
                organization="Humboldt-University Berlin, Institute for Theoretical Biology",
                site="https://livermetabolism.com",
            ),
            Creator(
                familyName="Höpfl",
                givenName="Sebastian",
                email="sebastian.hoepfl@isa.uni-stuttgart.de",
                organization="University of Stuttgart, Institute for Stochastics and Applications",
                site="https://github.com/shoepfl",
            ),
        ],
        units=U,
        model_units=ModelUnits(
            time=U.hr,
            extent=U.mmole,
            substance=U.mmole,
            length=U.meter,
            area=U.m2,
            volume=U.liter,
            # concentration [µg/l] -> [ng/ml]
        ),
        annotations=[
            # mouse (mus musculus)
            (BQB.IS, "snomed/447612001"),
            (BQB.IS, "ncit/C45247"),
            (BQB.IS, "taxonomy/10090"),
        ]
    )

    _m.compartments = [
        Compartment(
            sid="Vgut",
            name="gut compartment",
            value=0.050 * bodyweight,
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
            notes="""
            small fraction of total volume
            """
        ),
        Compartment(
            sid="Vcent",
            name="central compartment",
            value=0.080 * bodyweight,  # 80 [µl/g] * bodyweight [kg] * 1000 [g/kg];
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
            annotations=annotations.compartments["plasma"],
            notes="""
            The approximate blood volume of a mouse is 77-80 μl/g.
            """
        ),
        Compartment(
            sid="Vperi",
            name="peripheral compartment",
            value=(1.0 - 0.08 - 0.05) * bodyweight,
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
        ),
    ]

    _m.species = [
        Species(
            sid=f"{sid}_gut",
            name=f"{name} gut",
            compartment="Vgut",
            initialAmount=0,  # dose_bw * bodyweight / Mr
            hasOnlySubstanceUnits=False,
            substanceUnit=U.mmole,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
            notes="""
            handled in amount not concentration
            """,
        ),
        Species(
            sid=f"{sid}_cent",
            name=f"{name} plasma",
            compartment="Vcent",
            initialConcentration=0.0,
            substanceUnit=U.mmole,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
        ),
        Species(
            sid=f"{sid}_peri",
            name=f"{name} peripheral",
            compartment="Vperi",
            initialConcentration=0.0,
            substanceUnit=U.mmole,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species[sid],
        ),
    ]

    _m.parameters = [
        Parameter(
            f"Mr_{sid}", Mr, unit=U.g_per_mole,
            name=f"molecular weight {name}",
            sboTerm=SBO.MOLECULAR_MASS,
        ),
        Parameter(
        "conc_conversion", 1000.0, unit=U.ng_per_mug,
            name="conversion factor [µg/ml] -> [ng/ml]",
            # sboTerm=SBO.SYSTEMS_DESCRIPTION_CONSTANT,
        ),
        Parameter(
            f"{sid}_plasma", NaN, unit=U.ng_per_ml,
            name=f"{name} plasma [ng/ml]",
            constant=False,
            annotations=annotations.species[sid],
            # sboTerm=SBO.SYSTEMS_DESCRIPTION_PARAMETER,
        )
    ]

    _m.rules = [
        AssignmentRule(
            f"{sid}_plasma", f"{sid}_cent * Mr_{sid} * conc_conversion"  # [mmole/l]*[g/mole]=[µg/ml] -> [ng/ml]
        ),
    ]

    _m.reactions = [
        Reaction(
            sid="APPLICATION",
            name=f"application {name}",
            equation=f"-> {sid}_gut",
            formula=f"(166.989 * {dose_bw} * {bodyweight} / {Mr}) * (1 - exp(-time/tau))*exp(-time/tau)",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
            notes="""
            [mmole/hr]
            absorption from gut
            """,
            pars=[
                Parameter(
                    sid="tau",
                    name="response time",
                    value=0.012,
                    unit=U.hr,
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
        ),
        Reaction(
            sid="ABSORPTION",
            name=f"absorption {name}",
            equation=f"{sid}_gut -> {sid}_cent",
            formula=f"kabs * {sid}_gut",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
            notes="""
            [mmole/hr]
            absorption from gut
            """,
            pars=[
                Parameter(
                    sid="kabs",
                    name="absorption rate",
                    value=1E-2,
                    unit=U.l_per_hr,
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
        ),
        Reaction(
            sid="CLEARANCE",
            name=f"clearance {name}",
            equation=f"{sid}_cent ->",
            formula=f"CL * {sid}_cent",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
            notes="""
            [mmole/hr]
            clearance from central compartment
            """,
            pars=[
                Parameter(
                    sid="CL",
                    name=f"clearance {name}",
                    value=1E-2,
                    unit=U.l_per_hr,
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
        ),
        Reaction(
            sid="R1",
            name="transport peripheral (R1)",
            equation=f"{sid}_cent -> {sid}_peri",
            formula=f"Q * {sid}_cent",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
            notes="""
            [mmole/hr]
            distribution in peripheral compartment
            """,
            pars=[
                Parameter(
                    sid="Q",
                    name="distribution Q",
                    value=1E-1,
                    unit=U.l_per_hr,
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
        ),
        Reaction(
            sid="R2",
            name="transport central (R2)",
            equation=f"{sid}_peri -> {sid}_cent",
            formula=f"Q * {sid}_peri",
            sboTerm=SBO.BIOCHEMICAL_REACTION,
            notes="""
            [mmole/hr]
            distribution from peripheral compartment
            """,
        ),
    ]

    return _m


if __name__ == "__main__":

    @dataclass
    class ModelInfo:
        """Class for generating PK models."""
        sid: str
        name: str
        Mr: float
        dose_bw: float

    for info in [
        ModelInfo(
            sid="caf",
            name="caffeine",
            Mr=194.19,  # [g/mole]
            dose_bw=5  # [mg/kg]
        ),
        ModelInfo(
            sid="mid",
            name="midazolam",
            Mr=325.8,  # [g/mole]
            dose_bw=2  # [mg/kg]
        ),
        ModelInfo(
            sid="cod",
            name="codeine",
            Mr=299.4,  # [g/mole]
            dose_bw=2  # [mg/kg]
        )
    ]:
        _m = create_model_from_info(sid=info.sid, name=info.name, Mr=info.Mr, dose_bw=info.dose_bw)
        create_pk_model(
            model=_m,
            models_dir=os.path.join(Path(__file__).parent, 'results/'),
            equations=True,
            visualize=True,
            delete_session=False
        )



"""Caffeine PK model."""
from pathlib import Path

from sbmlutils.console import console
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.examples.templates import terms_of_use
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pk_models.models.templates import create_pk_model

# FIXME: [mmole/l] -> [ng/ml]
# FIXME: annotations of model components


class U(Units):
    """UnitDefinitions."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    mg = UnitDefinition("mg")
    per_min = UnitDefinition("per_min", "1/min")
    per_min_l = UnitDefinition("per_min_l", "1/min/liter")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/liter")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    l_per_min_mmole = UnitDefinition("l_per_min_mmole", "l/min/mmole")
    mM = UnitDefinition("mM", "mmole/liter")


_m = Model(
    sid="caffeine_pk",
    name="Caffeine pharmacokinetics model",
    notes="""
    # Model of absorption and distribution of caffeine.
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
    ],
    units=U,
    model_units=ModelUnits(
        time=U.min,
        extent=U.mmole,
        substance=U.mmole,
        length=U.meter,
        area=U.m2,
        volume=U.liter,
    ),
)
_m.compartments = [
    Compartment(
        sid="Vgut",
        name="gut compartment",
        value=1.0,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
    ),
    Compartment(
        sid="Vperi",
        name="peripheral compartment",
        value=1.0,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
    ),
    Compartment(
        sid="Vcent",
        name="central compartment",
        value=1.0,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
    ),
]

_m.species = [
    Species(
        sid="caf_gut",
        name="caffeine gut",
        compartment="Vgut",
        initialAmount=1.0,
        hasOnlySubstanceUnits=False,
        substanceUnit=U.mmole,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        notes="""
        handled in amount not concentration
        """,
    ),
    Species(
        sid="caf_cent",
        name="caffeine plasma",
        compartment="Vcent",
        initialConcentration=0.0,
        substanceUnit=U.mmole,
        sboTerm=SBO.SIMPLE_CHEMICAL,
    ),
    Species(
        sid="caf_peri",
        name="caffeine peripheral",
        compartment="Vperi",
        initialConcentration=0.0,
        substanceUnit=U.mmole,
        sboTerm=SBO.SIMPLE_CHEMICAL,
    ),
]

_m.reactions = [
    Reaction(
        sid="ABSORPTION",
        name="absorption caffeine",
        equation="caf_gut -> caf_cent",
        formula="k * caf_gut",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""
        [mmole/min]
        absorption from gut
        """,
        pars=[
            Parameter(
                sid="k",
                name="absorption rate",
                value=1.0,
                unit=U.l_per_min,
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
    ),
    Reaction(
        sid="CLEARANCE",
        name="clearance caffeine",
        equation="caf_cent ->",
        formula="CL * caf_cent",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""
        [mmole/min]
        clearance from central compartment
        """,
        pars=[
            Parameter(
                sid="CL",
                name="clearance caffeine",
                value=1.0,
                unit=U.l_per_min,
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
    ),
    Reaction(
        sid="R1",
        name="transport peripheral (R1)",
        equation="caf_cent -> caf_peri",
        formula="Q * caf_cent",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""
        [mmole/min]
        distribution in peripheral compartment
        """,
        pars=[
            Parameter(
                sid="Q",
                name="distribution Q",
                value=1.0,
                unit=U.l_per_min,
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
    ),
    Reaction(
        sid="R2",
        name="transport central (R2)",
        equation="caf_peri -> caf_cent",
        formula="Q * caf_peri",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        notes="""
        [mmole/min]
        distribution from peripheral compartment
        """,
    ),
]




if __name__ == "__main__":
    from pk_models import MODELS_DIR

    create_pk_model(model=_m, models_dir=MODELS_DIR, visualize=True)



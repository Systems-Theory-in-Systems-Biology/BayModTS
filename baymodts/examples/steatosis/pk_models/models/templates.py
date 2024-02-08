"""Reusable templates."""
from pathlib import Path

from sbmlutils.console import console
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import Creator, FactoryResult, Model, create_model
import os

creators = [
        Creator(
            familyName="König",
            givenName="Matthias",
            email="koenigmx@hu-berlin.de",
            organization="Humboldt-University Berlin, Institute for Theoretical Biology",
            site="https://livermetabolism.com",
        ),
    ]


def create_pk_model(model: Model, models_dir: Path, equations: bool = False, visualize: bool = False, delete_session: bool = True) -> None:
    """Create model."""
    results: FactoryResult = create_model(
        model=model,
        filepath=os.path.join(models_dir, f"{model.sid}.xml"),
        sbml_level=3,
        sbml_version=2,
    )

    # create differential equations
    if equations:
        md_path = os.path.join(models_dir, f"{model.sid}.md")
        ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
        ode_factory.to_markdown(md_file=md_path)

        console.rule(style="white")
        from rich.markdown import Markdown

        with open(md_path, "r") as f:
            md_str = f.read()
            md = Markdown(md_str)
            console.print(md)
        console.rule(style="white")

    # visualize network
    if visualize:
        visualize_sbml(sbml_path=results.sbml_path, delete_session=delete_session)

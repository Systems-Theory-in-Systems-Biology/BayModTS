from sbmlutils.metadata import *


compartments = {
    "plasma": [

    ]
}

species = {
    "caf": [
        # caffeine
        (BQB.IS, "pubchem.compound/2519"),
        (BQB.IS, "chebi/CHEBI:27732"),
        (BQB.IS, "ncit/C328"),
        (BQB.IS, "inchikey/RYYVLZVUVIJVGH-UHFFFAOYSA-N")
    ],
    "mid": [
        # midazolam
        (BQB.IS, "pubchem.compound/4192"),
        (BQB.IS, "chebi/CHEBI:6931"),
        (BQB.IS, "ncit/C62049"),
        (BQB.IS, "inchikey/DDLIGBOFAVUZHB-UHFFFAOYSA-N"),
        (BQB.IS, "snomed/373476007")
    ],
    "cod": [
        # codeine
        (BQB.IS, "pubchem.compound/5284371"),
        (BQB.IS, "chebi/CHEBI:16714"),
        (BQB.IS, "ncit/C383"),
        (BQB.IS, "inchikey/OROGSEYTTFOCAN-DNJOTXNNSA-N"),
        (BQB.IS, "snomed/387494007")
    ],
}
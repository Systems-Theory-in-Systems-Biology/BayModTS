# PK models for caffeine, midazolam and codeine
This folder contains the code to generate the PK models for the test substances caffeine, midazolam and codeine.

The drug cocktail consisted of caffeine, midazolam, and codeine (2 mg/kg body weight (bwt), 5 mg/kg bwt, and 2 mg/kg bwt, respectively). All in-vivo studies were carried out using male C57BL6/J mice (Janvier, France) 28–30 g of body
weight, between eight to ten months of age (ex-breeder) (n = 4–6/group). 

With body weight 29 g the following doses were applied:
- caffeine: 2 mg/kg * 0.029 kg = 0.058 mg
- midazolam: 5 mg/kg * 0.029 kg = 0.145 mg
- codeine: 2 mg/kg * 0.029 kg = 0.058 mg

## Installation
Create a virtual environment

```bash
mkvirtualenv baymodts --python=python3.11
cd baymodts/examples/steatosis/pk_models
pip install -r requirements.txt
```

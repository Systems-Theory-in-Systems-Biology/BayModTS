# model: codeine_pk
Autogenerated ODE System from SBML file with sbmlutils.
```
time: [hr]
substance: [mmol]
extent: [mmol]
volume: [l]
area: [m^2]
length: [m]
```

## Parameters `p`
```
CL = 0.01  # [l/hr] 
Mr_cod = 325.8  # [g/mol] 
Q = 0.1  # [l/hr] 
Vcent = 0.00232  # [l] 
Vgut = 0.00145  # [l] 
Vperi = 0.02523  # [l] 
conc_conversion = 1000.0  # [ng/µg] 
kabs = 0.01  # [l/hr] 
```

## Initial conditions `x0`
```
cod_cent = 0.0  # [mmol/l] Vcent
cod_gut = 0.12277470841006759  # [mmol/l] Vgut
cod_peri = 0.0  # [mmol/l] Vperi
```

## ODE system
```
# y
ABSORPTION = kabs * cod_gut  # [mmol/hr]
CLEARANCE = CL * cod_cent  # [mmol/hr]
R1 = Q * cod_cent  # [mmol/hr]
R2 = Q * cod_peri  # [mmol/hr]
cod_plasma = cod_cent * Mr_cod * conc_conversion  # [mg/ml]

# odes
d cod_cent/dt = (ABSORPTION / Vcent - CLEARANCE / Vcent - R1 / Vcent) + R2 / Vcent  # [mmol/l/hr]
d cod_gut/dt = -ABSORPTION / Vgut  # [mmol/l/hr]
d cod_peri/dt = R1 / Vperi - R2 / Vperi  # [mmol/l/hr]
```
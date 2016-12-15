# SGWB model selection project

## Installation

1. clone this repo
```python
cd sgwb_model_selection
pip install .
```

## Run a single power law model
```python
>>> from sbms.models import power_law
>>> omega_gw, f = power_law.omega_gw_spectrum(1e-9, alpha=3)
```

## Run a few power law models from a config file
Create a config file with specific parameters for your model (and make sure you use the right names...right now you'll have to look in the code for these). We can add a few power law models.

```ini
[power law]
omega_alpha=1e-9
alpha=3

[power law]
omega_alpha=1e-7
alpha=2/3
```

Run a model from a file:
```python
>>> from sbms.models import omgwf
>>> omega_gw, f = omgwf('my_ini_file.ini')
```

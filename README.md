# SGWB model selection project

[![Build Status](https://travis-ci.org/pmeyers279/sgwb_model_selection.svg?branch=master)](https://travis-ci.org/pmeyers279/sgwb_model_selection)

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
Create a config file with specific parameters for your model (and make sure you use the right names...right now you'll have to look in the code for these). We can add a few power law models. Right now fractions aren't supported for spectral indices. Sorry!

```ini
[power law_1]
omega_alpha=1e-9
alpha=3

[power law_2]
omega_alpha=1e-7
alpha=0.6667
```

*note* Right now the code checks if the second to last character is an underscore. If so it assumes that the model is specified by everything up to the underscore. This way we don't end up overwriting keys when we read in params.

Run a model from a file:
```python
>>> from sbms.models import omgwf
>>> omega_gw, f = omgwf('my_ini_file.ini')
```

# SGWB model selection project

[![Build Status](https://travis-ci.org/pmeyers279/sgwb_model_selection.svg?branch=master)](https://travis-ci.org/pmeyers279/sgwb_model_selection)

## Run a single power law model

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

## Installation

### CIT LIGO cluster

* Install PyMultiNest. See instructions [here](https://ldas-jobs.ligo.caltech.edu/~meyers/resources/resources/python/pymultinest.html)
* Install this package into the virtualenv by cloning the repo, cd to top level directory and then run `pip install .` or `pythons setup.py install`.
* It should...just work...

### OS X

1. Just don't do it right now?

2. If you need to...

#### Install this package

1. clone or fork this repo
```python
cd sgwb_model_selection
pip install .
```

#### Install MultiNest

2. Install gcc, cmake if not already installed. 

On OS X with macports:
```bash
sudo port install cmake gcc6
export FC=/opt/local/bin/gfortran-mp-6
```

3. Install MultiNest
```
git clone https://github.com/JohannesBuchner/MultiNest.git
cd MultiNest/build
cmake .. && make
```

4. Make sure that the libraries are available. Set these environment variables somewhere
```bash
export LAPACK=/usr/lib64/liblapack.so
export ATLAS=/usr/lib64/libatlas.so
export BLAS=/usr/lib64/libblas.so
export LD_LIBRARY_PATH=${MULTINEST_BASE}/MultiNest/lib/:$LD_LIBRARY_PATH
```
where `${MULTINEST_BASE}` is the base directory for wherever you downloaded multinest.

5. If you're on a mac...
You have to move a library to the lib/ of your python installation (and...change its name?), then set an environment variable properly. For me that was
```
mv ${MULTINEST_BASE}/lib/libmultinest.3.10.dylib ~/opt/sbms/lib/libmultinest.dylib
export DYLD_LIBRARY_PATH=~/opt/sbms/lib/
```

#### Install lalsimulation

This is needed for sensitivity curves right now (hopefully not for much longer...I was lazy).

Using macports: `sudo port install lalsimulation py27-lalsimulation


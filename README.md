# SGWB model selection project

[![Build Status](https://travis-ci.org/pmeyers279/sgwb_model_selection.svg?branch=master)](https://travis-ci.org/pmeyers279/sgwb_model_selection)

## Quick tutorial: running the code

Right now you just need to install the code in a virtualenv and then run `sbms-run [options]`. The possible flags are:

```bash
sbms-run -h
Usage: sbms-run [options]

Options:
  -h, --help            show this help message and exit
  -i INJ_FILE, --injection-file=INJ_FILE
                        injection param file
  -r RECOVERY_FILE, --recovery-file=RECOVERY_FILE
                        recovery param file
  -o OUTPUT_PREF, --output-prefix=OUTPUT_PREF
                        output prefix. path needs to exist right now.
```

### Injection file `-i`

The inj file specifies an injection to make. For right now you can't supply your data from a separate file (that should change relatively soon), but instead you just input a model for the stochastic background and the noise. There should always be at least a `[noise]` section with an observation time in years and a `noise_str`. (I'd recommend just using the one below in whatever files you want for now) and one model. If you want to inject multiple models of the same name, add an underscore with a number. An example of a file injecting a single model is below.

The header for the model should end in an `_#` and the rest of the header should match one of the available models (see below for current available models and their parameters).

```ini
[noise]
noise_str=-A
Tobs=1

[broken power law_1]
omg_ref=1e-6
alpha1=3
alpha2=-3
f_break=50
```

### Recovery File `-r`
This file is more confusing. Headers are just labels. There is a `model_type` parameter, which should be the name of a model in the registry with underscores replacing spaces. There is an `nparams` value taht should specify the number of variable parameters. Then there are up to 9 available slots for params. These params should be specified **in the order they are supplied to the model!!**. They are passed from the order they are found in the file straight into the model in python in the same order. See below for function calls for proper ordering.

Each param slot is a comma separated list. The first is the name that will appear in the plots and outputs from `PyMultiNest`. The second is a prior specification. See below for the key and available priors. The rest are inputs to the different priors.

**on the way** ability to specify that one parameter is always larger than another. This is nice when searching for multiple power laws or a broken power law of a specific form. I'm fairly certain I know how to implement this.

```
[1]
model_type=broken_power_law
nparams=4
;; 'name','prior type','input 1', 'input 2'
param1=omg_ref,LOG,1e-8,1e-5
param2=f_break,U,20,100
param3=alpha1,U,0,5
param4=alpha2,U,-5,0
```

### Available models:
Available models: 

Right now everything just runs from 20-100 Hz and for the power law fref defaults to 20 Hz. As soon as I work out how best to properly handle passing keyword arguments and optional arguments to the models, I'll add in the ability to change the frequency range.

Recall that the arguments below are explicit parameters set in the injection file and then for the recovery file, they are (in the order presented below) the params that need to be specified.

* `power law`
```
omg_alpha : `float`
    :math:`\Omega_{\alpha}`, amplitude for power spectrum
alpha : `float`
    spectral index
```
* `broken power law`
```
omg_ref : `float`
    :math:`\Omega_{\alpha}`, amplitude for power spectrum at the cross-over
    from one power law to the other
f_break : `float`
    frequency at which cross-over occurs
alpha1 : `float`
    power law index for first section
alpha2 : `float`
    power law index for second section
```



### Available Priors:
Type : `Key` [params]

(`Key` is the key for the recovery file)

* Uniform : `U` [start, end]
* Log : `LOG` [start, end]
* Delta : `DELTA` [return val, JUNK]
* Gaussian : `GAUSS` [mean, standard deviation]

### Example files
There should be examples of these files in the `examples/` directory.

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


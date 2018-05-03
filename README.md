# SLD to Volume Fractions Calculator

**Authors:** Filip Ciesielski, Luke Clifton from Rutherford Appleton Laboratory, STFC, Oxfordshire.


## Introduction

SLD to Volume Fractions Calculator is a simple script for inferring volume fractions of a single layer of mixed components.
At the moment only tripple-component layers are supported, by default comprising of a solvent, protein and lipid.

As of now, the library provides a function `find_volumes` that can be used in a custom python script, or in a jupyter notebook.
Web interface is planned for a near future as well, though.


## Installation Instructions

Requires python3

### Using virtualenv:

Clone the project folder into a local directory, then inside that folder:

```bash

$ virtualenv --python=python3 env
$ source ./env/bin/activate
(env)$ pip install -r requirements.txt
```

Then to launch the jupyter notebook:

```bash

(env)$ jupyter notebook
```

### Using pipenv:

Recommended way to install - using pipenv:

```bash

$ pipenv --three
$ pipenv install --dev
```

## Quickstart

The core of this library is based on the `find_volumes` function provided in the `helpers.core` module. Its use is 
demonstrated in a demo notebook called `example_notebook.ipynb`.
 
To launch the jupyter notebook:

### Launching notebook via virtualenv:

```bash
(env) > jupyter notebook 
```


### Launching notebook via pipenv: 
```bash
$ pipenv run jupyter notebook 
```

or 

```bash
$ pipenv shell            # <--- spawns the virtual env shell
$ jupyter notebook
```

### Using your own sample data:

User can modify the example notebook (it is a good idea to save a copy as a template) and provide your own values for

 
* `sample_d2o_conc`   (D2O fraction; 0<=, >=100)
* `sample_slds`     (SLD of the sample layer)
* `sample_sld_errors`  (SLD uncertainties)

Note, that all three have to have an equal number of values. Once filled, rerun the whole notebook.

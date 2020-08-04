![Build Status](https://travis-ci.com/JEstabrook/Enrich.svg?token=ZRDWBWe9sXCivP1NrZwq&branch=master)  [![](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-367) ![t](https://img.shields.io/badge/license-MIT-nrightgreen.svg) ![t](https://img.shields.io/badge/status-stable-nrightgreen.svg)


# Enrich


**regulon-enrichment** is a Python module used to predict the activity of regulatory proteins from RNAseq data.

Enrich submodules:

### `enricher.features` ###
Load -omic datasets


### `enricher.regulon` ###
Regulon utilities

# Dependencies

**regulon-enrichment** requires:
~~~~~~~~~~~~
- Python (>= 3.6)
- scikit-learn (>= 0.21.3)
- NumPy (>= 1.17.3)
- SciPy (>= 1.3.1)
- pandas (>= 0.25.3)
- tqdm (>= 4.38.0)
- dill (>= 0.3.1.1)
~~~~~~~~~~~~

# User installation
~~~~~~~~~~~~~~~~~

If you already have a working installation of numpy and scipy,
the easiest way to install regulon-enrichment is using ``pip``   ::

    pip install pip install regulon-enrichment==0.0.2b0

or ``conda``::

    conda install -c conda-forge regulon-enrichment
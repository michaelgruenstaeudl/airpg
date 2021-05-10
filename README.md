*airpg*: Automatically accessing the inverted repeats of archived plastid genomes
=================================================================================

[![Build Status](https://travis-ci.com/michaelgruenstaeudl/airpg.svg?branch=master)](https://travis-ci.com/michaelgruenstaeudl/airpg)
[![PyPI status](https://img.shields.io/pypi/status/airpg.svg)](https://pypi.python.org/pypi/airpg/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/airpg.svg)](https://pypi.python.org/pypi/airpg/)
[![PyPI version shields.io](https://img.shields.io/pypi/v/airpg.svg)](https://pypi.python.org/pypi/airpg/)
[![PyPI license](https://img.shields.io/pypi/l/airpg.svg)](https://pypi.python.org/pypi/airpg/)

A Python package for automatically accessing the inverted repeats of thousands of plastid genomes stored on NCBI Nucleotide

## INSTALLATION
To get the most recent stable version of *airpg*, run:

    pip install airpg

Or, alternatively, if you want to get the latest development version of *airpg*, run:

    pip install git+https://github.com/michaelgruenstaeudl/airpg.git


## EXAMPLE USAGE


### [Tutorial 1](airpg/tutorials/tutorial1.md): Very short survey (runtime ca. 5 min.; for the impatient)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide within the past 10 days.


### [Tutorial 2](airpg/tutorials/tutorial2.md): Short survey (runtime ca. 15 min.; for testing)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide within the current month.


### [Tutorial 3](airpg/tutorials/tutorial3.md): Medium survey (runtime ca. 5 hours)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide in 2019 only. Note: The results of this survey are available on Zenodo via DOI [10.5281/zenodo.4335906](https://zenodo.org/record/4335906)


### [Tutorial 4](airpg/tutorials/tutorial4.md): Full survey (runtime ca. 19 hours; with explanations)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide from January 2000 until, and including December 2020. Note: The results of this survey are available on Zenodo via DOI [10.5281/zenodo.4335906](https://zenodo.org/record/4335906)



<!--
## PACKAGING INSTRUCTIONS
```
#pip install .  ## For local testing

python3 -m build
python3 -m twine upload --repository testpypi dist/*
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps airpg

python3 -m twine upload dist/*
python3 -m pip install airpg
```
-->

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

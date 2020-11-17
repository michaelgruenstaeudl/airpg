*airpg*: Accessing the inverted repeats of archived plastid genomes
===================================================================

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
#### STEP 1: Generating plastome availability table
```
# Angiosperms
TESTFOLDER=./03_testing/angiosperms_Start2000toEnd2019
DATE=$(date '+%Y_%m_%d')
MYQUERY='complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 2000/01/01:2019/12/31[PDAT] AND 0000050000:00000250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND (Embryophyta[ORGN] AND Magnoliophyta[ORGN])'
AVAILTABLE=plastome_availability_table_${DATE}.tsv
mkdir -p $TESTFOLDER
```
```
# Defining blacklist
if [ ! -f ./BLACKLIST__master_${DATE} ]; then
    cat $(ls ./BLACKLIST__* | grep -v "master") > ./BLACKLIST__master_${DATE}
fi
```
```
airpg_retrieve.py -q "$MYQUERY" -o $TESTFOLDER/$AVAILTABLE --blacklist ./BLACKLIST__master_${DATE} 1>>$TESTFOLDER/airpg_retrieve_${DATE}.runlog 2>&1
```

#### STEP 2: Downloading records and extracting IR information
```
REPRTDSTAT=reported_IR_stats_table_${DATE}.tsv
mkdir -p $TESTFOLDER/records_${DATE}
mkdir -p $TESTFOLDER/data_${DATE}
```
```
airpg_analyze.py -i $TESTFOLDER/$AVAILTABLE -r $TESTFOLDER/records_${DATE}/ -d $TESTFOLDER/data_${DATE}/ -o $TESTFOLDER/$REPRTDSTAT 1>>$TESTFOLDER/airpg_analyze_${DATE}.runlog 2>&1
```

<!--
## FOO BAR BAZ
```
Foo bar baz
```
-->

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

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
#### STEP 1: Querying NCBI Nucleotide for complete plastid genomes given an Entrez search string
```
# Angiosperms
TESTFOLDER=./03_testing/angiosperms_Start2000toEnd2019
DATE=$(date '+%Y_%m_%d')
ENTREZSTRING='complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 2000/01/01:2019/12/31[PDAT] AND 0000050000:00000250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND (Embryophyta[ORGN] AND Magnoliophyta[ORGN])'
RECORDSTABLE=plastome_availability_table_${DATE}.tsv
mkdir -p $TESTFOLDER
```
```
# Updating blocklist
if [ ! -f ./airpg_blocklist.txt ]; then
    touch ./airpg_blocklist.txt
fi
airpg_update_blocklist.py -f ./airpg_blocklist.txt
```
```
airpg_identify.py -q "$ENTREZSTRING" -o $TESTFOLDER/$RECORDSTABLE \
    --blocklist ./airpg_blocklist.txt 1>>$TESTFOLDER/airpg_identify_${DATE}.runlog 2>&1
```

#### STEP 2: Retrieving and parsing the genome records identified in step 1, analyzing the position and length of their IR annotations
```
IRSTATSTABLE=reported_IR_stats_table_${DATE}.tsv
mkdir -p $TESTFOLDER/records_${DATE}
mkdir -p $TESTFOLDER/data_${DATE}
```
```
airpg_analyze.py -i $TESTFOLDER/$RECORDSTABLE \
    -r $TESTFOLDER/records_${DATE}/ -d $TESTFOLDER/data_${DATE}/ \
    -m john.smith@example.com -o $TESTFOLDER/$IRSTATSTABLE 1>>$TESTFOLDER/airpg_analyze_${DATE}.runlog 2>&1
```

<!--
## FOO BAR BAZ
```
Foo bar baz
```
-->

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

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

### Short survey for the impatient / for testing (runtime ca. 4 hours)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide in 2019 only. Note: The results of this survey are available on Zenodo via DOI [10.5281/zenodo.4335906](https://zenodo.org/record/4335906)
```
airpg_update_blocklist.py -f airpg_blocklist.txt \
-m john.smith@example.com -q "inverted[TITLE] AND \
repeat[TITLE] AND loss[TITLE]"

airpg_identify.py -q "complete genome[TITLE] AND \
(chloroplast[TITLE] OR plastid[TITLE]) AND \
2019/01/01:2019/12/31[PDAT] AND 50000:250000[SLEN] \
NOT unverified[TITLE] NOT partial[TITLE] AND \
(Embryophyta[ORGN] AND Magnoliophyta[ORGN])" \
-b airpg_blocklist.txt -o output_script1.tsv

airpg_analyze.py -i output_script1.tsv \
-m john.smith@example.com -o output_script2.tsv
```


### Full survey with explanations (runtime ca. 18 hours)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide from start of 2000 until end of October 2020. Note: The results of this survey are available on Zenodo via DOI [10.5281/zenodo.4335906](https://zenodo.org/record/4335906)

#### STEP 1: Querying NCBI Nucleotide for complete plastid genomes given an Entrez search string
```
TESTFOLDER=./angiosperms_Start2000toEndOct2020
DATE=$(date '+%Y_%m_%d')
ENTREZSTRING='complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 2000/01/01:2020/10/31[PDAT] AND 50000:250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND (Embryophyta[ORGN] AND Magnoliophyta[ORGN])' # complete plastid genomes of all flowering plants between start of 2000 and end of October 2020
RECORDSTABLE=plastome_availability_table_${DATE}.tsv
mkdir -p $TESTFOLDER
```
```
# Updating blocklist
if [ ! -f ./airpg_blocklist.txt ]; then
    touch ./airpg_blocklist.txt
    airpg_update_blocklist.py -f ./airpg_blocklist.txt
fi
airpg_update_blocklist.py -f ./airpg_blocklist.txt -m john.smith@example.com -q "inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]"
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

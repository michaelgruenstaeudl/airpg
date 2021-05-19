TUTORIAL 4
==========

#### Full survey (runtime ca. 19 hours; with explanations)

Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide from January 2000 until, and including, December 2020.

*Note*: The results of this survey are available on Zenodo via DOI [10.5281/zenodo.4335906](https://zenodo.org/record/4335906)


##### STEP 1: Querying NCBI Nucleotide for complete plastid genomes given an Entrez search string
```
TESTFOLDER=./angiosperms_Start2000toEndOct2020
DATE=$(date '+%Y_%m_%d')
ENTREZSTRING='complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 2000/01/01:2020/12/31[PDAT] AND 50000:250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND Magnoliopsida[ORGN]' # complete plastid genomes of all flowering plants between start of 2000 and end of 2020
RECORDSTABLE=plastome_availability_table_${DATE}.tsv
mkdir -p $TESTFOLDER

# Generating (or updating) blocklist
airpg_update_blocklist.py -f ./airpg_blocklist.txt -m john.smith@example.com -q "inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]"

airpg_identify.py -q "$ENTREZSTRING" -o $TESTFOLDER/$RECORDSTABLE \
    --blocklist ./airpg_blocklist.txt 1>>$TESTFOLDER/airpg_identify_${DATE}.runlog 2>&1

# Sorting output_script1.tsv and output_script1.tsv.duplicates by UID number
sort -k1 -n output_script1.tsv > tmp && mv tmp output_script1.tsv
sort -k1 -n output_script1.tsv.duplicates > tmp && mv tmp output_script1.tsv.duplicates
```

##### STEP 2: Retrieving and parsing the genome records identified in step 1, analyzing the position and length of their IR annotations
```
IRSTATSTABLE=reported_IR_stats_table_${DATE}.tsv
mkdir -p $TESTFOLDER/records_${DATE}
mkdir -p $TESTFOLDER/data_${DATE}

airpg_analyze.py -i $TESTFOLDER/$RECORDSTABLE \
    -r $TESTFOLDER/records_${DATE}/ -d $TESTFOLDER/data_${DATE}/ \
    -m john.smith@example.com -o $TESTFOLDER/$IRSTATSTABLE 1>>$TESTFOLDER/airpg_analyze_${DATE}.runlog 2>&1
```

##### STEP 3: Confirming the presence, position, and length of the IR annotations identified in step 2 through self-blasting of each sequence record via [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/).
```
EXTENDEDIRSTATS=extended_IR_stats_table_${DATE}.tsv
mkdir -p $TESTFOLDER/records_${DATE}
mkdir -p $TESTFOLDER/data_${DATE}

airpg_confirm.py -i $TESTFOLDER/$IRSTATSTABLE \
    -d $TESTFOLDER/data_${DATE}/ -n 10000 -x 40000 \
    -o $TESTFOLDER/$EXTENDEDIRSTATS 1>>$TESTFOLDER/airpg_confirm_${DATE}.runlog 2>&1
```

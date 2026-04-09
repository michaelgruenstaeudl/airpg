#!/bin/bash


#TUTORIAL 1

# Very short survey (runtime ca. 5 min.; for the impatient)

# Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide within the past 10 days.

# Place '-k api_key' in arguments if you have one.

STARTDATE=$(date -d "10 days ago" +%Y/%m/%d)
ENDDATE=$(date +%Y/%m/%d)

airpg_identify.py \
-q "complete genome[TITLE] AND \
(chloroplast[TITLE] OR plastid[TITLE]) AND \
$STARTDATE:$ENDDATE[PDAT] AND \
50000:250000[SLEN] NOT unverified[TITLE] \
NOT partial[TITLE] AND Magnoliopsida[ORGN]" \
-o output_script1.tsv \
-e john.smith@example.com
#&> output_script1.log

mkdir -p records
mkdir -p data

airpg_analyze.py \
-i output_script1.tsv \
-m john.smith@example.com \
-o output_script2.tsv \
--recordsdir records/ \
--datadir data/ \
#&> output_script2.log

airpg_confirm.py \
-i output_script2.tsv \
-o output_script3.tsv \
-n 10000 -x 50000 \
--datadir data/ \
#&> output_script3.log
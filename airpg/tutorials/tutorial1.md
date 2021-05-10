TUTORIAL 1
==========

##### Very short survey (runtime ca. 5 min.; for the impatient)

Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide within the past 10 days.


```
TODAY=$(date +%d)
if (($TODAY >= 6 && $TODAY <= 10)); then
    STARTDATE=$(date +%Y/%m/01)
elif (($TODAY >= 11 && $TODAY <= 15)); then
    STARTDATE=$(date +%Y/%m/05)
elif (($TODAY >= 16 && $TODAY <= 20)); then
    STARTDATE=$(date +%Y/%m/10)
elif (($TODAY >= 21 && $TODAY <= 25)); then
    STARTDATE=$(date +%Y/%m/15)
else
    PREVMONTH=$(printf "%02d" $(($(date +%m)-1)))
    STARTDATE=$(date +%Y/$PREVMONTH/20)
fi
ENDDATE=$(date +%Y/%m/%d)

airpg_identify.py \
-q "complete genome[TITLE] AND \
(chloroplast[TITLE] OR plastid[TITLE]) AND \
$STARTDATE:$ENDDATE[PDAT] AND \
50000:250000[SLEN] NOT unverified[TITLE] \
NOT partial[TITLE] AND Magnoliopsida[ORGN]" \
-o output_script1.tsv \
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
```

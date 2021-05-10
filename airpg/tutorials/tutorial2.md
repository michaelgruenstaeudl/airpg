TUTORIAL 2
==========
Short survey (runtime ca. 15 min.; for testing)
Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide within the current month.

```
airpg_identify.py -q "complete genome[TITLE] AND \
(chloroplast[TITLE] OR plastid[TITLE]) AND \
$(date +%Y/%m/01):$(date +%Y/%m/%d)[PDAT] AND \
50000:250000[SLEN] NOT unverified[TITLE] \
NOT partial[TITLE] AND Magnoliopsida[ORGN]" \
-o output_script1.tsv # &> output_script1.log

airpg_analyze.py -i output_script1.tsv \
-m john.smith@example.com -o output_script2.tsv \
# &> output_script2.log
```

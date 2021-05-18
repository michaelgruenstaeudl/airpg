TUTORIAL 3
==========

#### Medium survey (runtime ca. 5 hours)

Survey of all plastid genomes of flowering plants submitted to NCBI Nucleotide in 2019 only.

*Note*: The results of this survey are available on Zenodo via DOI [10.5281/zenodo.4335906](https://zenodo.org/record/4335906)


```
airpg_update_blocklist.py -f airpg_blocklist.txt \
-m john.smith@example.com -q "inverted[TITLE] AND \
repeat[TITLE] AND loss[TITLE]"

airpg_identify.py -q "complete genome[TITLE] AND \
(chloroplast[TITLE] OR plastid[TITLE]) AND \
2019/01/01:2019/12/31[PDAT] AND 50000:250000[SLEN] \
NOT unverified[TITLE] NOT partial[TITLE] AND \
Magnoliopsida[ORGN]" \
-b airpg_blocklist.txt -o output_script1.tsv

airpg_analyze.py -i output_script1.tsv \
-m john.smith@example.com -o output_script2.tsv
```

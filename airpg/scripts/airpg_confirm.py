#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a complete plastid genome sequence (in FASTA format) and re-calculates (re-infers) the position (and, thus, the length) of the IRs. The input is the output of script 02. The output is a table of plastid genome records (one record per row) that lists the originally inferred IR length and the BLAST-inferred IR length so that a comparison is possible.

DESIGN:
    * Like in the other scripts, the evaluation of the genome records is conducted one by one, not all simultaneously.

    * The following batch code shall be used to re-calculate the compare the IRs (i.e., the IR file pair) of each record:
        ```
        # Self-blasting of the plastid genome sequence in order to infer the IR length
        blastn -db chloroplastGenome.fasta -query chloroplastGenome.fasta -outfmt 7 -strand 'both' | awk '{ if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}'
        ```

TO DO (FUTURE VERSIONS OF AIRPG):
    * Start positions for IRa and IRb are sometimes switched around (probably because they were improperly recorded in previous scripts due to missing unique identification)

    * Once the IRs are re-calculated, the originally inferred IR length and the newly calculated IR length shall be compared to see if previous studies have - on average - overestimated or underestimated the IR length.

    * If differences between the originally inferred IR length and the newly calculated IR length are discovered, it will be interesting to see on which side of the IRs (the side that faces the LSC or the side that faces the SSC) the original inference was incorrect (i.e., on which side a bias in the original inference happened).

    * Also do the inference of the IR via MUMMER (specifically, a self-comparison via dnadiff) so that the IR boundaries as inferred via self-BLASTING are confirmed (i.e., similar to the internal confirmation check of the total number of sequence differences via CMP).

'''

#####################
# IMPORT OPERATIONS #
#####################
import os, argparse
import subprocess
import pandas as pd
import coloredlogs, logging

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019-2021 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Retrieve the plastid genomes identified by the first script '\
           'and evaluate their inverted repeats'
__version__ = '2021.05.18'

#############
# DEBUGGING #
#############
#import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############

## THE FOLLOWING LINES CAN BE IMPLEMENTED IN A FUTURE VERSION OF AIRPG:
#def coerceToExactLocation(location):
#    exactLocation = None
#    if '<' in str(location) or '>' in str(location):
#        exactLocation = str(location)[1:]
#    else:
#        exactLocation = str(location)
#    return exactLocation


def main(args):
  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    if args.verbose:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)
    else:
        coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='INFO', logger=log)

  # STEP 2. Read table from Script02 and list accessions
    try:
        IR_table = pd.read_csv(args.infn, index_col=False, sep='\t', encoding="utf-8", na_values="n.a.")
        accessions = list(IR_table["ACCESSION"].values)
        IR_table = IR_table.astype({
            "ACCESSION": str,
            "IRa_REPORTED": str,
            "IRa_REPORTED_START": pd.Int64Dtype(),
            "IRa_REPORTED_END": pd.Int64Dtype(),
            "IRa_REPORTED_LENGTH": pd.Int64Dtype(),
            "IRb_REPORTED": str,
            "IRb_REPORTED_START": pd.Int64Dtype(),
            "IRb_REPORTED_END": pd.Int64Dtype(),
            "IRb_REPORTED_LENGTH": pd.Int64Dtype()
        })
    except Exception as err:
        log.exception("Error reading file `%s`: %s" % (str(args.infn), str(err)))
        raise

  # STEP 3. Add new columns to table for accession
    added_columns = [
        "IRa_BLASTINFERRED",
        "IRa_BLASTINFERRED_START",
        "IRa_BLASTINFERRED_END",
        "IRa_BLASTINFERRED_LENGTH",
        "IRb_BLASTINFERRED",
        "IRb_BLASTINFERRED_START",
        "IRb_BLASTINFERRED_END",
        "IRb_BLASTINFERRED_LENGTH"]

    if not any(col in list(IR_table.columns) for col in added_columns):
        IR_table = IR_table.reindex(columns = list(IR_table.columns) + added_columns)

    IR_table = IR_table.astype({
        "IRa_BLASTINFERRED": str,
        "IRa_BLASTINFERRED_START": pd.Int64Dtype(),
        "IRa_BLASTINFERRED_END": pd.Int64Dtype(),
        "IRa_BLASTINFERRED_LENGTH": pd.Int64Dtype(),
        "IRb_BLASTINFERRED": str,
        "IRb_BLASTINFERRED_START": pd.Int64Dtype(),
        "IRb_BLASTINFERRED_END": pd.Int64Dtype(),
        "IRb_BLASTINFERRED_LENGTH": pd.Int64Dtype()
        })

    IR_table = IR_table.set_index("ACCESSION", drop=True)

  # STEP 4. Loop over accession in inlist
    # Step 4.1. Check if FASTA file for accession exists
    main_dir = os.getcwd()
    for accession in accessions:
        try:
            acc_folder = os.path.join(args.datadir, str(accession))
            seq_FASTA = accession + "_completeSeq.fasta"
        except Exception as err:
            log.warning("Error accessing FASTA file of accession `%s`: %s.\nSkipping this accession." % (str(accession), str(err)))
            continue

        # Step 4.2. Change into accession folder and conduct BLAST locally
        # Change to directory containing sequence files
        os.chdir(acc_folder)
        # Set up local BLAST database
        filestem_db = accession + "_completeSeq" + "_blastdb"
        try:
            log.info("Creating local BLAST database for accession `%s`." % (str(accession)))
            mkblastargs = ["makeblastdb", "-in", seq_FASTA, "-parse_seqids", "-title", accession, "-dbtype", "nucl", "-out", filestem_db, "-logfile", filestem_db+".log"]
            mkblastdb_subp = subprocess.Popen(mkblastargs)
            returncode = mkblastdb_subp.wait()
            if returncode != 0: # Can probably be done prettier
                raise Exception
        except Exception as err:
            log.exception("Error creating local BLAST database for accession `%s`: %s\nSkipping this accession." % (str(accession), str(err)))
            continue

        # Infer IR positions through self-BLASTing
        try:
            log.info("Self-BLASTing FASTA file of accession `%s` to identify the IRs." % (str(accession)))
            blastargs = ["blastn", "-db", filestem_db, "-query", seq_FASTA, "-outfmt", "7", "-strand", "both"]
            blast_subp = subprocess.Popen(blastargs, stdout=subprocess.PIPE)
            awkargs = ["awk", "{if ($4 > " + str(args.minlength) + " && $4 < " + str(args.maxlength) + ") print $4, $7, $8, $9, $10}"]
            awk_subp = subprocess.Popen(awkargs, stdin=blast_subp.stdout, stdout=subprocess.PIPE)
            out, err = awk_subp.communicate()
            result_lines = out.splitlines()
        except Exception as err:
            log.warning("Error while self-BLASTing FASTA file of accession `%s`: %s.\nSkipping this accession." % (str(accession), str(err)))
            continue

        # Compress local BLAST database if BLAST output received
        if len(result_lines) != 0:
            try:
                tarargs = ["tar", "czf", filestem_db+"_FILES.tar.gz", filestem_db+".*", "--remove-files"] # "--remove-files" must be at end
                returncode = subprocess.call(" ".join(tarargs), shell=True) # Shell=True is necessary for the wildcard
                if returncode != 0:                                     # Can probably be done prettier
                    raise Exception("Non-zero exit status")    # TO BE IMPROVED: # Error message of subprocess.call is not transferred to exception
            except Exception as err:
                log.warning("Error while compressing local BLAST database for accession `%s`: %s." % (str(accession), str(err)))

        # Step 4.3. Parse output of self-BLASTing
        # Note: BLAST sometimes finds additional regions in the sequence that match the length requirements filtered for in awk. We only want the IRs, and therefore need to pick out the two regions with matching length
        if len(result_lines) > 2:
            temp_lines = []
            for i in range(len(result_lines)-1):
                for j in range(i+1,len(result_lines)):
                    if result_lines[i].split()[1] == result_lines[j].split()[1]:
                        temp_lines.append(result_lines[i])
                        temp_lines.append(result_lines[j])
                        break
            result_lines = temp_lines

        if len(result_lines) == 2:
            # Compare the start positions of the found regions. By default, we assume IRb is located before IRa in the sequence
            if result_lines[0].split()[1] > result_lines[1].split()[1]:
                IRa_info = result_lines[1].split()
                IRb_info = result_lines[0].split()
            else:
                IRa_info = result_lines[0].split()
                IRb_info = result_lines[1].split()

        # Step 4.4. Save data into correct columns
        # Note: It is important to stay within the condition 'len(result_lines) == 2'
            
            IR_table.at[accession, "IRa_BLASTINFERRED"] = "yes"
            IR_table.at[accession, "IRb_BLASTINFERRED"] = "yes"

            IR_table.at[accession, "IRa_BLASTINFERRED_START"] = int(IRa_info[1])
            IR_table.at[accession, "IRb_BLASTINFERRED_START"] = int(IRb_info[1])

            IR_table.at[accession, "IRa_BLASTINFERRED_END"] = int(IRa_info[2])
            IR_table.at[accession, "IRb_BLASTINFERRED_END"] = int(IRb_info[2])

            IR_table.at[accession, "IRa_BLASTINFERRED_LENGTH"] = int(IRa_info[0])
            IR_table.at[accession, "IRb_BLASTINFERRED_LENGTH"] = int(IRb_info[0])

            ## THE FOLLOWING LINES CAN BE IMPLEMENTED IN A FUTURE VERSION OF AIRPG:
            #            IR_table.at[accession, "IRa_START_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRa_REPORTED_START"])) - float(coerceToExactLocation(IR_table.at[accession, "IRa_BLASTINFERRED_START"])))
            #            IR_table.at[accession, "IRb_START_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRb_REPORTED_START"])) - float(coerceToExactLocation(IR_table.at[accession, "IRb_BLASTINFERRED_START"])))

            #            IR_table.at[accession, "IRa_END_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRa_REPORTED_END"])) - float(coerceToExactLocation(IR_table.at[accession, "IRa_BLASTINFERRED_END"])))
            #            IR_table.at[accession, "IRb_END_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRb_REPORTED_END"])) - float(coerceToExactLocation(IR_table.at[accession, "IRb_BLASTINFERRED_END"])))

            #            IR_table.at[accession, "IRa_LENGTH_COMPARED_DIFFERENCE"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRa_REPORTED_LENGTH"])) - float(coerceToExactLocation(IR_table.at[accession, "IRa_BLASTINFERRED_LENGTH"])))
            #            IR_table.at[accession, "IRb_LENGTH_COMPARED_DIFFERENCE"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRb_REPORTED_LENGTH"])) - float(coerceToExactLocation(IR_table.at[accession, "IRb_BLASTINFERRED_LENGTH"])))

        else:
            log.info("Could not infer IRs for accession `%s`:\n%s." % (str(accession), "\n".join([str(line).strip() for line in result_lines])))    # TO BE IMPROVED: Does not always work!
            IR_table.at[accession, "IRa_BLASTINFERRED"] = "no"
            IR_table.at[accession, "IRb_BLASTINFERRED"] = "no"

        # Step 4.5. Change back to main directory
        os.chdir(main_dir)

        # TO BE IMPROVED: It would be best if the script would append to the output file after each accession, and not write the full output file after processing all accessions.

  # STEP 5. Save extended IR list to outfile
    IR_table.to_csv(args.outfn, sep='\t', header=True, na_rep="n.a.")


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--infn", "-i", type=str, required=True, help="Path to input file; input is a summary table on reported IR positions and length (tab-delimited, accession numbers in first column)")
    parser.add_argument("--outfn", "-o", type=str, required=True, help="Path to output file that contains extended table IR positions and length")
    parser.add_argument("--datadir", "-d", type=str, required=True, default="./data/", help="Path to folder containing record-specific subfolders that store each record's complete sequence in FASTA format")
    parser.add_argument("--minlength", "-n", type=int, required=False, default="10000", help="(Optional) Minimal length of IR for BLASTing")
    parser.add_argument("--maxlength", "-x", type=int, required=False, default="50000", help="(Optional) Maximum length of IR for BLASTing")
    parser.add_argument("--verbose", "-v", action="store_true", required=False, default=False, help="(Optional) Enable verbose logging")
    args = parser.parse_args()
    main(args)

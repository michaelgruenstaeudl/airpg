#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a complete plastid genome sequence (in FASTA format) and re-calculates (re-infers) the position (and, thus, the length) of the IRs.

    The output shall be a table of plastid genome records (one record per row) that lists the originally inferred IR length and this newly BLASTINFERRED IR length so that a comparison is possible.

DESIGN:
    * Like in the other scripts, the evaluation of the records is conducted one by one, not all simultaneously.

    * The plastid genome sequence of each record will be bundled with the GB file of that record in a record-specific gzip file. Hence, before that record is evaluated here, that specific gzip-file must be unpacked; after the evaluation, the gzip-file must be reconstituted.

TO DO (FOR REVISION):
    * Start positions for IRa and IRb are sometimes switched around (probably because they were improperly recorded in previous scripts due to missing unique identification)

TO DO (LONG-TERM):
    * Once the IRs are re-BLASTINFERRED, the originally inferred IR length and the newly BLASTINFERRED IR length shall be compared in order to see if previous studies have - on average - overestimated or underestimated the IR length.

    * If differences between the originally inferred IR length and the newly BLASTINFERRED IR length are discovered, it will be interesting to see on which side of the IRs (the side that faces the LSC or the side that faces the SSC) the original inference was incorrect (i.e., on which side a bias in the original inference happened).

    * The following batch code shall be used to re-calculate the compare the IRs (i.e., the IR file pair) of each record:
        ```
        # Self-blasting of the plastid genome sequence in order to infer the IR length
        blastn -db chloroplastGenome.fasta -query chloroplastGenome.fasta -outfmt 7 -strand 'both' | awk '{ if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}'
        ```

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
__version__ = '2021.05.10'

#############
# DEBUGGING #
#############
#import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############

def coerceToExactLocation(location):
    exactLocation = None
    if '<' in str(location) or '>' in str(location):
        exactLocation = str(location)[1:]
    else:
        exactLocation = str(location)
    return exactLocation

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
        "IRb_BLASTINFERRED_LENGTH": pd.Int64Dtype(),
        })

  # STEP 4. Loop over accession in inlist
    # Step 4.1. Check if FASTA file for accession exists
    main_dir = os.getcwd()
    for accession in accessions:
        try:
            acc_folder = os.path.join(args.datadir, str(accession))
            seq_FASTA = os.path.join(acc_folder, accession + "_completeSeq.fasta")
            # Test if file exists: os.path.isfile(os.path.join(acc_folder, accession + "_completeSeq.fasta"))
        except Exception as err:
            log.warning("Error accessing FASTA file of accession `%s`: %s.\nSkipping this accession." % (str(accession), str(err)))
            continue

        # Step 4.2. Change into accession folder and conduct BLAST locally
        # Change to directory containing sequence files
        os.chdir(acc_folder)
        # Set up local BLAST database
        try:
            mkblastargs = ["makeblastdb", "-in", seq_FASTA, "-parse_seqids", "-title", accession, "-dbtype", "nucl"]
            mkblastdb_subp = subprocess.Popen(mkblastargs)
            mkblastdb_subp.wait()
        except Exception as err:
            log.exception("Error creating blastdb for accession `%s`: %s\nSkipping this accession." % (str(accession), str(err)))
            continue
        # Infer IR positions through self-BLASTing
        try:
            blastargs = ["blastn", "-db", seq_FASTA, "-query", seq_FASTA, "-outfmt", "7", "-strand", "both"]
            blast_subp = subprocess.Popen(blastargs, stdout=subprocess.PIPE)
            awkargs = ["awk", "{if ($4 >", str(args.minlength), "&& $4 <", str(args.minlength), ") print $4, $7, $8, $9, $10}"]
            awk_subp = subprocess.Popen(awkargs, stdin=blast_subp.stdout, stdout=subprocess.PIPE)
            out, err = awk_subp.communicate()
            result_lines = out.splitlines()
            print(result_lines)
        except Exception as err:
            log.warning("Error while self-BLASTing FASTA file of accession `%s`: %s.\nSkipping this accession." % (str(accession), str(err)))
            continue

'''
    # append columns that will be filled in this script
    # Step 2.2. Get accessions



    #main_dir = os.getcwd()
    #for folder in folders:
    #    # Init values that will be written to table
    #    accession = os.path.basename(folder)


        try:
            blastargs = ["blastn", "-db", str(accession) + "_completeSeq.fasta", "-query", accession + "_completeSeq.fasta", "-outfmt", "7", "-strand", "both"]
            blast_subp = subprocess.Popen(blastargs, stdout=subprocess.PIPE)
            awkargs = ["awk", "{if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}"]
            awk_subp = subprocess.Popen(awkargs, stdin=blast_subp.stdout, stdout=subprocess.PIPE)
            out, err = awk_subp.communicate()
            result_lines = out.splitlines()
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
                    ira_info = result_lines[1].split()
                    irb_info = result_lines[0].split()
                else:
                    ira_info = result_lines[0].split()
                    irb_info = result_lines[1].split()

                # Assign BLASTINFERRED info
                IR_table.at[accession, "IRa_BLASTINFERRED"] = "yes"
                IR_table.at[accession, "IRb_BLASTINFERRED"] = "yes"

                IR_table.at[accession, "IRa_BLASTINFERRED_START"] = int(ira_info[1])
                IR_table.at[accession, "IRb_BLASTINFERRED_START"] = int(irb_info[1])

                IR_table.at[accession, "IRa_BLASTINFERRED_END"] = int(ira_info[2])
                IR_table.at[accession, "IRb_BLASTINFERRED_END"] = int(irb_info[2])

                IR_table.at[accession, "IRa_BLASTINFERRED_LENGTH"] = int(ira_info[0])
                IR_table.at[accession, "IRb_BLASTINFERRED_LENGTH"] = int(irb_info[0])

                IR_table.at[accession, "IRa_START_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRa_REPORTED_START"])) - float(coerceToExactLocation(IR_table.at[accession, "IRa_BLASTINFERRED_START"])))
                IR_table.at[accession, "IRb_START_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRb_REPORTED_START"])) - float(coerceToExactLocation(IR_table.at[accession, "IRb_BLASTINFERRED_START"])))

                IR_table.at[accession, "IRa_END_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRa_REPORTED_END"])) - float(coerceToExactLocation(IR_table.at[accession, "IRa_BLASTINFERRED_END"])))
                IR_table.at[accession, "IRb_END_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRb_REPORTED_END"])) - float(coerceToExactLocation(IR_table.at[accession, "IRb_BLASTINFERRED_END"])))

                IR_table.at[accession, "IRa_LENGTH_COMPARED_DIFFERENCE"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRa_REPORTED_LENGTH"])) - float(coerceToExactLocation(IR_table.at[accession, "IRa_BLASTINFERRED_LENGTH"])))
                IR_table.at[accession, "IRb_LENGTH_COMPARED_DIFFERENCE"] = int(float(coerceToExactLocation(IR_table.at[accession, "IRb_REPORTED_LENGTH"])) - float(coerceToExactLocation(IR_table.at[accession, "IRb_BLASTINFERRED_LENGTH"])))

            else:
                log.warning("Could not calculate IRs for accession " + accession + "." + "\n".join([str(line).strip() for line in result_lines]))

                IR_table.at[accession, "IRa_BLASTINFERRED"] = "no"
                IR_table.at[accession, "IRb_BLASTINFERRED"] = "no"

        except Exception as err:
            log.exception("Error while calculating IRs. %s\n Skipping this accession." % (str(err)))
            continue

'''

        # WHY SHOULD THIS BE APPENDED EVERY ACCESSION AND NOT AT THE END?
        # Step 4.x. Write data to outfn
        IR_table.to_csv(args.outfn, sep='\t', header=True, na_rep="n.a.")
        os.chdir(main_dir)


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--infn", "-i", type=str, required=True, help="Path to input file; input is a summary table on reported IR positions and length (tab-delimited, accession numbers in first column)")
    parser.add_argument("--outfn", "-o", type=str, required=True, help="Path to output file that contains extended table IR positions and length")
    parser.add_argument("--datadir", "-d", type=str, required=True, help="Data folder containing subfolders with FASTA files")
    parser.add_argument("--minlength", "-n", type=int, required=False, default="10000", help="(Optional) Minimal length of IR for BLASTing")
    parser.add_argument("--maxlength", "-x", type=int, required=False, default="50000", help="(Optional) Maximum length of IR for BLASTing")
    parser.add_argument("--verbose", "-v", action="store_true", required=False, default=False, help="(Optional) Enable verbose logging")
    args = parser.parse_args()
    main(args)

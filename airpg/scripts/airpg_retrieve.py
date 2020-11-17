#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
	This script generates a table of plastid genome records that are currently available on NCBI. It simultaneously ensures that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records). If an output file already exists, this script will identify the most recent date of any of the processed records and search NCBI only for younger records.

	The output is a table in which each row contains the parsed information of a single record. Each row contains nine, tab-separated columns in the following order:
	01. the unique identifier,
	02a. the accession number,
	02b. synonyms of the accession number (i.e., issue about regular vs. RefSeq NC_ records),
	03. the sequence version number,
	04. the organism name,
	05. the sequence length,
	06. the date the record went online,
	07. the authors (uppermost AUTHORS line in GB-file),
	08. the name of the publication (uppermost TITLE line in GB-file), and
	09. the full citation of the publication (see uppermost JOURNAL line in GB-file)
	10. Any note if a REFSEQ accession number exists and to which regular accession number it is equal to.
	11. the taxonomy

DESIGN:
	There are thousands of plastid genome sequences on GenBank. The parsing of the records is, thus, conducted one by one, not all simultaneously. Specifically, a list of unique identifiers is first obtained and then this list is looped over.

TO DO:
	* none for now

NOTES:
	* none for now

'''

#####################
# IMPORT OPERATIONS #
#####################
import os.path
import argparse
import coloredlogs, logging
#import pandas as pd
from airpg import entrez_interaction
from airpg import table_io
from datetime import datetime

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
			 'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019-2020 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Collect summary information on all plastid sequences stored ' \
		   'in NCBI Nucleotide'
__version__ = '2020.08.31.1930'

#############
# DEBUGGING #
#############
#import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############


def main(args):

  # STEP 1. Set up logger
	log = logging.getLogger(__name__)
	coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

	EI = entrez_interaction.EntrezInteraction(log)

  # STEP 2. Check if output file already exists, read existing UIDs, infer mindate
	uids_already_processed = []
	min_date = None
	outfn = os.path.abspath(args.outfn)

	if args.blacklist:
		tio = table_io.TableIO(outfn, fp_blacklist = args.blacklist, logger = log)
	else:
		tio = table_io.TableIO(outfn, logger = log)
	fp_duplicates = os.path.join(os.path.dirname(outfn), os.path.basename(outfn) + ".duplicates")
	if os.path.isfile(fp_duplicates):
		tio.read_duplicates(fp_duplicates)
		uids_already_processed.extend(map(int, tio.duplicates.keys()))

	if len(tio.entry_table) > 0:
		uids_already_processed.extend(list(tio.entry_table.index.values))
		log.info("Summary file '%s' already exists. Number of UIDs read: %s" % (str(outfn), str(len(uids_already_processed))))
		if args.update_only:
			min_date = datetime.strptime(tio.entry_table["CREATE_DATE"].max(), '%Y-%m-%d')
			log.info("NOTE: Only records more recent than '%s' are being looked for." % (str(min_date)))
	else:
		log.info(("Summary file '%s' does not exist; generating new file. Thus, no UIDs read." % (str(outfn))))

	# STEP 3. Get all existing UIDs and calculate which to be processed
	try:
		uids_new = EI.retrieve_uids(args.query, min_date)
	except Exception as err:
		log.exception("Error while retrieving UID list: " + str(err))
	log.info(("Number of UIDs on NCBI: %s" % (str(len(uids_new)))))
	uids_to_process = set(uids_new) - set(uids_already_processed)
	log.info(("Number of UIDs to be processed: %s" % (str(len(uids_to_process)))))

	# STEP 4. Parse all entries, append entry-wise to file
	for uid in uids_to_process:
		log.info(("Reading and parsing UID '%s', writing to '%s'." % (str(uid), str(outfn))))
		try:
			xml_entry = EI.fetch_xml_entry(uid)
			parsed_entry = EI.parse_xml_entry(xml_entry)
		except Exception as err:
			log.exception("Error retrieving info for UID " + str(uid) + ": " + str(err) + "\nSkipping this accession.")
			continue
		duplseq = parsed_entry.pop("DUPLSEQ")
		tio.entry_table.loc[uid] = parsed_entry
		if duplseq:
			tio.duplicates[uid] = [parsed_entry["ACCESSION"], duplseq]
		tio.append_entry_to_table(parsed_entry, uid, outfn)

	# STEP 5. Remove duplicates of REFSEQs and blacklisted entries
	tio.write_duplicates(fp_duplicates)
	tio.remove_duplicates()
	tio.remove_blacklisted_entries()
	'''
	# Alternative way (i.e. with logs) to remove duplicates
	for refseq, dup in tio.duplicates.items():
		try:
			tio.entry_table.drop(tio.entry_table.loc[tio.entry_table["ACCESSION"] == dup].index, inplace = True)
			log.info("Removed accession '%s' from '%s' because it is a duplicate of REFSEQ '%s'." % (dup, str(os.path.basename(outfn)), refseq))
		except: # If the duplicate doesn't exist, nothing happens (except for a log message)
			log.info("Could not find accession '%s' when trying to remove it." % str(dup))
	'''
	tio.write_entry_table(outfn)

########
# MAIN #
########

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
	parser.add_argument("-o", "--outfn", type=str, required=True, help="path to output file")
	parser.add_argument("-q", "--query", type=str, required=False, default="complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 0000050000:00000250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND (Embryophyta[ORGN] AND Magnoliophyta[ORGN])", help="(Optional) Entrez query that will replace the standard query")
	parser.add_argument("-u", "--update_only", action="store_true", required=False, default=False, help="Will only add entries with more recent creation date than the most recent existing entry.")
	parser.add_argument("-b", "--blacklist", type=str, required=False, help="path to file of blacklisted genera that will be removed from the retrieved plastid sequences.")
	args = parser.parse_args()
	main(args)

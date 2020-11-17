#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
	This script takes a list of GenBank accession numbers and - for each accession number - downloads the record
	from GenBank, parses it via Biopython, extracts the inverted repeat (IR) regions and saves both their reported
	position as well as their reported sequences to files. The IR is hereby identified either explicitly (i.e., via
	features with appropriate keywords) or implicitly (i.e., via annotations of the large and small single
	copy region and then taking the complement set thereof).

	The output shall be a set of files for each GenBank record:
	(a) the GenBank record in GB format
	(b) the full plastid genome sequence in FASTA format
	(c) the reported positions of IRa and IRb as a table
	(d) the reported sequence of IRa and IRb (reverse-complemented), both in FASTA format

	The downloaded GB record is gzipped and saved in a folder titled "records".
	The other files generated (i.e., the three FASTA files and the table) are saved in a folder called "data" and,
	within that, their own subfolder (named with the accession number)

DESIGN:
	* The IRs (i.e. IRa and IRb) of a plastid genome record may be annotated in different ways (i.e., via different
	  features and feature qualifiers), depending on the record. This script is flexible enough to identify the
	  different naming conventions, yet always extract only a single IR pair per record.

	* This script can also infer the position of the IR implicitly via the position of large (LSC) and small single
	  copy (SSC) regions. Specifically, it extracts the positions of the LSC and the SSC and calculate the IRs as
	  the complement set of the LSC and the SSC, if no explicit annotation information for the IR is present.

	* The plastid genome records on GenBank are inconsistent in their annotations of the reading direction of the
	  IR. Specifically, the annotations sometimes indicate that one of the IRs is on the opposite strand than its
	  counterpart (i.e., is reverse-complemented compared to its counterpart), sometimes not. In order to avoid
	  issues when extracting the IRs in FASTA-format, an explicit check of the IRs per genome, followed by a
	  reverse-complementing operation (where required), is indicated. The situation is complicated by the fact
	  that a certain level of non-identity (i.e., up to 10% of all nucleotides) must remain permissible between
	  the IRs, because the identification of such a non-identity is the overall objective of this investigation.
	  Consequently, this script explicitly checks for each GenBank record (in function "getInvertedRepeats()") if
	  one of the IR features of a record must be reverse complemented or not. Specifically, the string conducts
	  approximate string comparisons and, based on the results, switches the strand information for one of the
	  IRs where necessary.

	* The output files of each record shall be bundled together as a record-specific gzip file.

TO DO:
	* see code locations indicated with "# TO DO"

NOTES:
	* none for now
'''

#####################
# IMPORT OPERATIONS #
#####################
from Bio import SeqIO
from fuzzywuzzy import fuzz
from ete3 import NCBITaxa
from pathlib import Path
from airpg import entrez_interaction
from airpg import table_io
from airpg import article_mining
from airpg import ir_operations
import pandas as pd
import os, argparse
import tarfile, coloredlogs, logging
import time

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
			 'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019-2020 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Compare IRs for a series of IR FASTA files'
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
	if args.verbose:
		coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)
	else:
		coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='INFO', logger=log)
	mail = args.mail
	query = args.query
	iro = ir_operations.IROperations(log)
	EI = entrez_interaction.EntrezInteraction(log)

  # STEP 2. Read in accession numbers to loop over
	tio = table_io.TableIO(args.infn, args.outfn, args.blacklist, logger = log)
	tio.remove_blacklisted_entries()

	accessions = list(tio.entry_table["ACCESSION"].values)
	if len(accessions) > 0:
		if not os.path.exists(args.recordsdir):
			os.makedirs(args.recordsdir)
		if not os.path.exists(args.datadir):
			os.makedirs(args.datadir)

  # STEP 3. Loop over accession in inlist
	for accession in accessions:
		acc_folder = os.path.join(args.datadir, str(accession))
		if not os.path.exists(acc_folder):
			os.makedirs(acc_folder)
		else:
			log.warning("Folder for accession `%s` already exists. Skipping this accession." % (str(accession)))
			continue

		log.info("Saving GenBank flat file for accession `%s`." % (str(accession)))
		try:
			fp_entry = EI.fetch_gb_entry(accession, acc_folder)
		except:
			log.warning("Error retrieving accession `%s`. Skipping this accession." % (str(accession)))
			os.rmdir(acc_folder)
			continue

		try:
			try:
				rec = SeqIO.read(fp_entry, "genbank")
			except Exception as err:
				raise Exception("Error reading record of accession `%s`: `%s`. Skipping this accession." %
				(str(accession), str(err)))
				continue

			rec_id = str(rec.id).split('.')[0]
			# Note: This internal check ensures that we are actually dealing with the record that was
			# intended to be downloaded via efetch.
			if not rec_id == str(accession):
				log.warning("Accession number mismatch. Expected: `%s`. Retrieved: `%s`. Skipping this accession." % \
				  (str(accession), rec_id))
				continue
			log.info("Writing sequence as FASTA for accession `%s`." % (str(accession)))
			iro.write_sequence_to_fasta(str(rec.seq), ">" + str(accession) + "_completeSequence", os.path.join(acc_folder, str(accession) + "_completeSeq.fasta"))

			ira_feature = None
			irb_feature = None
			if not str(accession) in tio.ir_table.index:
				tio.ir_table = tio.ir_table.append(pd.Series(name=str(accession)))
			try:
				ira_feature, irb_feature = iro.identify_inverted_repeats(rec, 1000)
				rev_comp = False
				if ira_feature and irb_feature:
					score_noRC = fuzz.ratio(ira_feature.extract(rec).seq,
											irb_feature.extract(rec).seq)
					score_RC = fuzz.ratio(ira_feature.extract(rec).seq,
										  irb_feature.extract(rec).seq.reverse_complement())
					if score_noRC < score_RC:
						rev_comp = True
				ir_info = iro.collect_info_from_features(ira_feature, irb_feature)
				tio.ir_table.loc[accession] = ir_info
				tio.append_ir_info_to_table(ir_info, accession, args.outfn)
			except Exception as err:
				ir_info = iro.collect_info_from_features(ira_feature, irb_feature)
				tio.ir_table.loc[accession] = ir_info
				tio.append_ir_info_to_table(ir_info, accession, args.outfn)
				raise Exception("Error while extracting IRs for accession `%s`: `%s`. Skipping further processing of this accession." % (str(accession), str(err)))
				continue
			tio.ir_table.loc[accession] = iro.collect_info_from_features(ira_feature, irb_feature)
			# TODO: Currently, nothing is done with the table object, since the file is written entry-wise. Remove full table from use?
			iro.write_irs_to_fasta(rec, ira_feature, irb_feature, acc_folder, rev_comp)
		except Exception as err:
			log.warning(str(err))
		finally:
			tar = tarfile.open(os.path.join(args.recordsdir, accession + ".tar.gz"), "w:gz")
			tar.add(fp_entry, os.path.basename(fp_entry))
			tar.close()
			os.remove(fp_entry)

  # STEP 4. Check any accession for IR loss and remove from outlist if necessary
	am = article_mining.ArticleMining(log)
	articles = EI.fetch_pubmed_articles(mail, query)
	ncbi = NCBITaxa()
	# Update database if it is older than 1 month
	if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 2592000:
		ncbi.update_taxonomy_database()
	article_genera = set()
	for article in articles:
		article_genera.union(am.get_genera_from_pubmed_article(article, ncbi))
	tio.read_ir_table(args.outfn)
	tio.remove_naturally_irl_genera(article_genera)
	tio.write_ir_table(args.outfn)


########
# MAIN #
########

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
	parser.add_argument("--infn", "-i", type=str, required=True, help="path to input file; input is a summary table of NCBI accessions (tab-delimited, accession numbers in second column)")
	parser.add_argument("--outfn", "-o", type=str, required=True, help="path to output file that contains information on IR positions and length")
	parser.add_argument("--recordsdir", "-r", type=str, required=False, default="./records/", help="path to records directory")
	parser.add_argument("--datadir", "-d", type=str, required=False, default="./data/", help="path to data directory")
	parser.add_argument("--verbose", "-v", action="store_true", required=False, default=False, help="Enable verbose logging.")
	parser.add_argument("--blacklist", "-b", type=str, required=False, help="path to taxonomy blacklist")
	parser.add_argument("--query", "-q", type=str, required=False, default="inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]", help="query to find pubmed articles describing inverted repeat loss")
	parser.add_argument("--mail", "-m", type=str, required=True, help="Mail account for PubMed entrez search. Any valid mail address will work.")
	args = parser.parse_args()
	#if bool(args.query) ^ bool(args.mail):
	#	parser.error("--query and --mail must be given together")
	main(args)

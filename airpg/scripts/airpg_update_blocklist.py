#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
	This script updates a taxon blocklist with genus names of taxa that have been indicated to lack inverted repeats in their plastid genomes through an automated query of NCBI PubMed. Interacts with NCBI Taxonomy to confirm the validity of genus names.

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
import time
#import PlastomeIntegrityChecks as pic
from ete3 import NCBITaxa
from pathlib import Path
from datetime import datetime

# For suppressing console output
import io
from contextlib import redirect_stdout

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019-2020 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Append a list of genus names of taxa that have been indicated to lack one or more inverted repeats in their plastid genome through an automated query of NCBI PubMed'
__version__ = '2020.12.17.1200'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

def get_irl_clade_species(ncbi):
    species_ids = []
    irl_clade_tree = ncbi.get_topology([ncbi.get_name_translator(['IRL clade'])['IRL clade'][0]])
    for leaf in irl_clade_tree.iter_leaves():
         species_ids.append(int(leaf.name))
    species = set(ncbi.translate_to_names(species_ids))
    return species

def get_irl_clade_genera(ncbi):
    genera_ids = []
    irl_clade_tree = ncbi.get_topology([ncbi.get_name_translator(['IRL clade'])['IRL clade'][0]])
    for node in irl_clade_tree.iter_descendants():
        if ncbi.get_rank([int(node.name)])[int(node.name)] == "genus":
            genera_ids.append(int(node.name))
    genera = set(ncbi.translate_to_names(genera_ids))
    return genera

def read_blocklist(fp_blocklist):
    '''
    Read a file of blocklisted genera/species.
    Params:
     - fp_blocklist: file path to input file
    '''
    blocklist = set()
    with open(fp_blocklist, "r") as fh_blocklist:
        for line in [l.strip() for l in fh_blocklist.readlines()]:
            if not line.startswith("#"):
                blocklist.add(line)
    return blocklist

def append_blocklist(fp_blocklist, blocklist):
    with open(fp_blocklist, "a") as fh_blocklist:
        fh_blocklist.write("# Update on %s\n" % datetime.now().strftime("%Y-%m-%d, %H:%M"))
        for entry in blocklist:
            fh_blocklist.write(entry + "\n")

def main(args):

    ## STEP 1. Initialize variables
    ncbi = NCBITaxa()
    # Update database if it is older than one month
    if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 2592000:
        ncbi.update_taxonomy_database()
    blocklist = set()
    blocklist_existing = set()

    ## STEP 2. Read blocklist if the file already exists
    if os.path.isfile(args.file_blocklist):
        print("Reading existing blocklist ...")
        blocklist_existing = read_blocklist(args.file_blocklist)

    ## STEP 3. Assemble species names of IRL clade of Fabaceae
    print("\nFetching genus names of taxa in 'IRL clade' of Fabaceae ...")
    irl_clade_genera = get_irl_clade_genera(ncbi)
    print("  Adding new species names to blocklist ...")
    blocklist = irl_clade_genera.difference(blocklist_existing)

    ## STEP 4. Append only new taxon names to blocklist
    print("\nCalculating and appending species names not previously in blocklist ...")
    append_blocklist(args.file_blocklist, blocklist)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-f", "--file_blocklist", type=str, required=True, help="path to blocklist file")
    parser.add_argument("-q", "--query", type=str, required=False, default="inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]", help="query used to fetch PMC articles that will be scanned for taxa with missing IRs")
    args = parser.parse_args()
    main(args)
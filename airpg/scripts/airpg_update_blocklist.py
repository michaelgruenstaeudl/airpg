#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script generates or updates a taxon blocklist with genus names of taxa that have been indicated to lack inverted repeats in their plastid genomes through an automated query of NCBI PubMed. Interacts with NCBI Taxonomy to confirm the validity of genus names.

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
from airpg import table_io
from airpg import article_mining as AM
from airpg import entrez_interaction as EI
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
__copyright__ = 'Copyright (C) 2019-2021 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Write a list of genus names of taxa that have been indicated to lack inverted repeats in their plastid genomes through an automated query of NCBI PubMed to the blocklist'
__version__ = '2021.05.19.1400'

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

def write_blocklist(fp_blocklist, blocklist):
    with open(fp_blocklist, "w") as fh_blocklist:
        fh_blocklist.write("# Blocklist as of %s\n" % datetime.now().strftime("%Y-%m-%d, %H:%M"))
        for entry in blocklist:
            fh_blocklist.write(entry + "\n")

def append_blocklist(fp_blocklist, blocklist):
    with open(fp_blocklist, "a") as fh_blocklist:
        fh_blocklist.write("# Update on %s\n" % datetime.now().strftime("%Y-%m-%d, %H:%M"))
        for entry in blocklist:
            fh_blocklist.write(entry + "\n")

def main(args):

    ## STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='INFO', logger=log)

    ## STEP 2. Initialize variables
    ncbi = NCBITaxa()
    # Update database if it is older than one month
    if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 2592000:
        ncbi.update_taxonomy_database()
    blocklist = set()
    blocklist_existing = set()

    ## STEP 3. Read blocklist if the file already exists
    if os.path.isfile(args.file_blocklist):
        log.info("Reading existing blocklist ...")
        blocklist_existing = read_blocklist(args.file_blocklist)

    ## STEP 4a. Collect genus names of IRL clade of Fabaceae
    log.info("Fetching genus names of taxa in 'IRL clade' of Fabaceae ...")
    try:
        irl_clade_genera = get_irl_clade_genera(ncbi)
    except:
        irl_clade_genera = set()    
    log.info("Adding new genus names to blocklist ...")
    blocklist = set(list(irl_clade_genera.union(blocklist_existing)))

    ## STEP 4b. Collect species names of IRL clade of Fabaceae
    log.info("Fetching species names of taxa in 'IRL clade' of Fabaceae ...")
    try:
        irl_clade_species = get_irl_clade_species(ncbi)
    except:
        irl_clade_species = set()    
    log.info("Adding new species names to blocklist ...")
    blocklist = set(list(irl_clade_species.union(blocklist)))
    
    ## STEP 5. Conduct the search on NCBI PubMed
    if args.query and args.mail:
        log.info("Querying NCBI Pubmed for taxon names ...")
        try:
            irl_clade_genera = set()
            am = AM.ArticleMining(log)
            ei = EI.EntrezInteraction(log)
            if ei.internet_on():  # Check if internet connection active
                articles = ei.fetch_pubmed_articles(args.mail, args.query)
            else:  # If no internet connection, raise error
                raise Exception("ERROR: No internet connection.")
            ncbi = NCBITaxa()
            # Update database if it is older than 1 month
            if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 2592000:
                ncbi.update_taxonomy_database()
            article_genera = set()
            for article in articles:
                article_genera.union(am.get_genera_from_pubmed_article(article, ncbi))
        except:
            article_genera = set()    
        blocklist = blocklist.union(article_genera)

    ## STEP 6. Keeping only species name if corresponding genus name present
    log.info("Removing genus names if individual species of genus in list ...")
    species_names = list()
    genus_names = set()
    genus_epithets = set()
    for line in blocklist:
        first, *rest = line.split()
        if rest:
            species_names.append(line)
            genus_epithets.add(first)
        else:
            genus_names.add(first)
    species_names.extend(genus_names - genus_epithets)
    blocklist = set(species_names)  # Making sure that no two lines are identical

    ## STEP 7. Write updated blocklist to file or replace old blocklist
    log.info("Writing updated blocklist to file ...")
    #append_blocklist(args.file_blocklist, blocklist)
    write_blocklist(args.file_blocklist, sorted(blocklist))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-f", "--file_blocklist", type=str, required=True, help="Path to blocklist file")
    parser.add_argument("-q", "--query", type=str, required=False, default="inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]", help="(Optional) Entrez string to query NCBI PubMed")
    parser.add_argument("-m", "--mail", type=str, required=False, help="(Optional) Mail address needed for Entrez search on NCBI PubMed (any valid mail address works)")
    args = parser.parse_args()
    main(args)

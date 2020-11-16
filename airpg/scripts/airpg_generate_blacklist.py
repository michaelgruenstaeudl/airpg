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
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Create or append a list of species names that are proven to lack one or more inverted repeats in their plastid genome'
__version__ = '2020.08.20.1800'

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

def read_blacklist(fp_blacklist):
    '''
    Read a file of blacklisted genera/species.
    Params:
     - fp_blacklist: file path to input file
    '''
    blacklist = set()
    with open(fp_blacklist, "r") as fh_blacklist:
        for line in [l.strip() for l in fh_blacklist.readlines()]:
            if not line.startswith("#"):
                blacklist.add(line)
    return blacklist

def append_blacklist(fp_blacklist, blacklist):
    with open(fp_blacklist, "a") as fh_blacklist:
        fh_blacklist.write("# Update on %s\n" % datetime.now().strftime("%Y-%m-%d, %H:%M"))
        for entry in blacklist:
            fh_blacklist.write(entry + "\n")

def main(args):

    ## STEP 1. Initialize variables
    ncbi = NCBITaxa()
    # Update database if it is older than one month
    if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 2592000:
        ncbi.update_taxonomy_database()
    blacklist = set()
    blacklist_existing = set()

    ## STEP 2. Read blacklist if the file already exists
    if os.path.isfile(args.file_blacklist):
        print("Reading existing blacklist ...")
        blacklist_existing = read_blacklist(args.file_blacklist)

    ## STEP 3. Assemble species names of IRL clade of Fabaceae
    print("\nFetching genus names of taxa in 'IRL clade' of Fabaceae ...")
    irl_clade_genera = get_irl_clade_genera(ncbi)
    print("  Adding new species names to blacklist ...")
    blacklist = irl_clade_genera.difference(blacklist_existing)

    ## STEP 4. Append only new taxon names to blacklist
    print("\nCalculating and appending species names not previously in blacklist ...")
    append_blacklist(args.file_blacklist, blacklist)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-f", "--file_blacklist", type=str, required=True, help="path to blacklist file")
    parser.add_argument("-q", "--query", type=str, required=False, default="inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]", help="query used to fetch PMC articles that will be scanned for species with missing IRs")
    args = parser.parse_args()
    main(args)

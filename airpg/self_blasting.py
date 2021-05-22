import shutil
import subprocess
import logging

class SelfBlasting:

    def __init__(self, seq_FASTA, accession, logger = None):
        if not shutil.which("blastn"):
            raise Exception("Error: 'blastn' not installed!")
        self.log = logger or logging.getLogger(__name__ + ".SelfBlasting")
        self.seq_FASTA = seq_FASTA
        self.accession = accession
        self.filestem_db = self.accession + "_completeSeq" + "_blastdb"        
        
    def setup_blast_db(self):
        mkblastargs = ["makeblastdb", "-in", self.seq_FASTA, "-parse_seqids", "-title", self.accession, "-dbtype", "nucl", "-out", self.filestem_db, "-logfile", self.filestem_db + ".log"]
        mkblastdb_subp = subprocess.Popen(mkblastargs)
        returncode = mkblastdb_subp.wait()
        if returncode != 0: # Can probably be done prettier
            raise Exception
        
    def infer_irs(self, minlength, maxlength):
        blastargs = ["blastn", "-db", self.filestem_db, "-query", self.seq_FASTA, "-outfmt", "7", "-strand", "both"]
        blast_subp = subprocess.Popen(blastargs, stdout=subprocess.PIPE)
        awkargs = ["awk", "{if ($4 > " + str(minlength) + " && $4 < " + str(maxlength) + ") print $4, $7, $8, $9, $10}"]
        awk_subp = subprocess.Popen(awkargs, stdin=blast_subp.stdout, stdout=subprocess.PIPE)
        out, err = awk_subp.communicate()
        return out.splitlines()
        
    def compress_db(self):
        tarargs = ["tar", "czf", self.filestem_db+"_FILES.tar.gz", self.filestem_db+".*", "--remove-files"] # "--remove-files" must be at end
        returncode = subprocess.call(" ".join(tarargs), shell=True) # Shell=True is necessary for the wildcard
        if returncode != 0:                                     # Can probably be done prettier
            raise Exception("Non-zero exit status")             # Error message of subprocess.call is not transferred to exception (because shell=True is set above, but the latter is necessary)

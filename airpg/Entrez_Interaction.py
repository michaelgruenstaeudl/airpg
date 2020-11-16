import os, subprocess, logging
import xml.etree.ElementTree as ET
from airpg.airpg import fetch_pubmed
import entrezpy.conduit
from ete3 import NCBITaxa
from datetime import date

class EntrezInteraction:

    def __init__(self, logger = None):
        self.log = logger or logging.getLogger(__name__ + ".EntrezInteraction")

    def retrieve_uids(self, query, min_date = None):
        '''
        Retrieves and returns a list of UIDs that are the result of an Entrez query
        Params:
         - query: the Entrez query
         - min_date: date. the minimum data to start searching from.
        '''
        if min_date:
            esearch_args = ['esearch', '-db', 'nucleotide', '-sort', '"Date Released"', '-mindate', min_date.strftime("%Y/%m/%d"), '-maxdate', date.today().strftime("%Y/%m/%d"), '-query', query]
        else:
            esearch_args = ['esearch', '-db', 'nucleotide', '-sort', '"Date Released"', '-query', query]
        esearch = subprocess.Popen(esearch_args, stdout=subprocess.PIPE)
        efetchargs = ["efetch", "-db", "nucleotide", "-format", "uid"]
        efetch = subprocess.Popen(efetchargs, stdin=esearch.stdout, stdout=subprocess.PIPE)
        out, err = efetch.communicate()

        return list(map(int, out.splitlines()))[::-1]

    def fetch_xml_entry(self, uid):
        '''
        Fetches one GenBank entry in GenBank XML format and returns it as ElementTree
        Params:
         - uid: the UID of the GenBank entry
        Returns: ElementTree of entry
        '''
        self.log.debug("Fetching XML GenBank entry " + str(uid))
        esummaryargs = ["esummary", "-db", "nucleotide", "-format", "gb", "-mode", "xml", "-id", str(uid)]
        esummary = subprocess.Popen(esummaryargs, stdout=subprocess.PIPE)
        out, err = esummary.communicate()
        return ET.fromstring(out)

    def parse_xml_entry(self, entry):
        '''
        Parses a GenBank XML entry and returns a dictionary of field name/value pairs
        Params:
         - entry: ElementTree. The XML entry
        '''
        # Parse out the relevant info from XML-formatted record summary
        uid_data = entry.find("GBSeq")
        fields = {}
        accession = uid_data.find("GBSeq_primary-accession").text
        fields["ACCESSION"] = accession
        fields["VERSION"] = uid_data.find("GBSeq_accession-version").text.split('.')[1]
        fields["ORGANISM"] = uid_data.find("GBSeq_organism").text
        fields["SEQ_LEN"] = uid_data.find("GBSeq_length").text
        fields["TAXONOMY"] = uid_data.find("GBSeq_taxonomy").text

        # Parse and format the date that the record was first online
        month_map = {"JAN":"01", "FEB":"02", "MAR":"03", "APR":"04", "MAY":"05", "JUN":"06", "JUL":"07", "AUG":"08", "SEP":"09", "OCT":"10", "NOV":"11", "DEC":"12"}
        create_date = uid_data.find("GBSeq_create-date").text.split('-')
        fields["CREATE_DATE"] = create_date[2] + "-" + month_map[create_date[1]] + "-" + create_date[0]

        # Parse all info related to the authors and the publication
        references = uid_data.find("GBSeq_references").findall("GBReference")
        authstring = ""
        title = ""
        citation = ""

        # Properly format the author output
        for ref in references:
            # Look for a reference that has authors (not all entries have a reference with authors)
            authors = ref.find("GBReference_authors")
            if authors:
                title = ref.find("GBReference_title").text
                citation = ref.find("GBReference_journal").text
                authors = ref.find("GBReference_authors").findall("GBAuthor")
                for author in authors:
                    authstring = authstring + author.text + ","
                authstring = authstring[:-2]
                break
        fields["AUTHORS"] = authstring
        fields["TITLE"] = title
        fields["REFERENCE"] = citation

        # Parse comment field for RefSeq
        note = ""
        duplseq = None
        if accession[:3] == "NC_":
            comments = uid_data.find("GBSeq_comment").text
            for comment in comments.split(";"):
                keyw1 = "PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review"
                keyw2 = "The reference sequence is identical to"
                keyw3 = "REVIEWED REFSEQ: This record has been curated by NCBI staff."
                keyw4 = "The reference sequence was derived from"
                keywords = [keyw1, keyw2, keyw3, keyw4]
                if any(keyw in comment for keyw in keywords): 
                    duplseq = comment.split(" ")[-1][:-1]
                    note = "The REFSEQ accession '%s' is identical to accession '%s'." % (accession, duplseq)
        fields["NOTE"] = note
        fields["DUPLSEQ"] = duplseq

        return fields

    def fetch_gb_entry(self, acc_id, outdir):
        '''
        Saves GenBank flatfile with accession number "acc_id" to
        outdir and returns the path and filename of the file
        Params:
         - acc_id: accession number of the GenBank entry
         - outdir: file path to output directory
        '''
        self.log.debug("Fetching GenBank entry %s and saving to %s" % (str(acc_id), outdir))
        gbFile = os.path.join(outdir, str(acc_id) + ".gb")
        with open(gbFile, "w") as outfile:
            efetchargs = ["efetch", "-db", "nucleotide", "-format", "gb", "-id", str(acc_id)]
            efetch = subprocess.Popen(efetchargs, stdout=outfile)
            efetch.wait()
        if not os.path.isfile(gbFile):
            raise Exception("Error retrieving GenBank flatfile of accession " + str(acc_id))
        elif os.path.getsize(gbFile) == 0:
            raise Exception("Error retrieving GenBank flatfile of accession " + str(acc_id))
        return gbFile

    def fetch_pubmed_articles(self, mail, query):
        '''
        Fetches all articles from PubMed found by query and returns them as a list of PubMedRecord objects
        Params:
         - mail: Mail address of requester
         - query: Entrez search string
        '''
        articles = None
        cond = entrezpy.conduit.Conduit(mail)
        fetch_pipe = cond.new_pipeline()
        sid = fetch_pipe.add_search({'db': 'pubmed', 'term': query, 'rettype': 'count'})
        fid = fetch_pipe.add_fetch({'retmode':'xml'}, dependency=sid, analyzer=fetch_pubmed.PubMedAnalyzer())
        a = cond.run(fetch_pipe)
        result = a.get_result()
        if result.size() >= 1:
            articles = result.pubmed_records
        return articles

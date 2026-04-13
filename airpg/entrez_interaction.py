import os, subprocess, logging, urllib.request, time
import xml.etree.ElementTree as ET
from airpg import parse_pubmed
import entrezpy.conduit
from ete3 import NCBITaxa
from datetime import date
from Bio import Entrez

class EntrezInteraction:

    def __init__(self, email, api_key=None, logger = None):
        self.email = email
        self.api_key = api_key
        self.log = logger or logging.getLogger(__name__ + ".EntrezInteraction")

    @property
    def _sleep_interval(self):
        return 0.11 if self.api_key else 0.4

    def retrieve_uids(self, query, min_date=None):
        '''
        Retrieves and returns a list of UIDs that are the result of an Entrez query
        Params:
         - query: the Entrez query
         - min_date: date. the minimum date to start searching from.
        '''
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key

        search_params = {
            "db": "nucleotide",
            "term": query,
            "sort": "Date Released",
            "usehistory": "y",
            "retmax": 100000
        }

        if min_date:
            search_params["mindate"] = min_date.strftime("%Y/%m/%d")
            search_params["maxdate"] = date.today().strftime("%Y/%m/%d")
            search_params["datetype"] = "pdat"

        handle = Entrez.esearch(**search_params)
        record = Entrez.read(handle)
        handle.close()

        return record["IdList"]

    def fetch_xml_entry(self, uid):
        '''
        Fetches one GenBank entry in GenBank XML format and returns it as ElementTree
        Params:
         - uid: the UID of the GenBank entry
        Returns: ElementTree of entry
        '''
        self.log.debug("Fetching XML GenBank entry " + str(uid))
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        handle = Entrez.efetch(db="nucleotide", id=str(uid), rettype="gb", retmode="xml")
        out = handle.read()
        handle.close()
        return ET.fromstring(out).find("GBSeq")  # Return GBSeq directly

    def parse_xml_entry(self, entry):
        '''
        Parses a GenBank XML entry and returns a dictionary of field name/value pairs
        Params:
         - entry: ElementTree. The XML entry
        '''
        # Parse out the relevant info from XML-formatted record summary
        uid_data = entry
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
        note = None
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
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        handle = Entrez.efetch(db="nucleotide", id=str(acc_id), rettype="gb", retmode="text")
        with open(gbFile, "w") as outfile:
            outfile.write(handle.read())
        handle.close()
        if not os.path.isfile(gbFile):
            raise Exception("Error retrieving GenBank flatfile of accession " + str(acc_id))
        elif os.path.getsize(gbFile) == 0:
            raise Exception("Error retrieving GenBank flatfile of accession " + str(acc_id))
        return gbFile

    def fetch_gb_entries_batch(self, acc_ids, outdir, batch_size=25):
        '''
        Fetches multiple GenBank flatfiles in batches and saves them to outdir.
        Returns a dict of {acc_id: filepath} for successfully downloaded entries.
        '''
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key

        results = {}
        acc_ids = list(acc_ids)
        MAX_RETRIES = 3
        RETRY_DELAY = 5  # seconds

        for i in range(0, len(acc_ids), batch_size):
            batch = acc_ids[i:i + batch_size]
            self.log.debug("Fetching batch %d-%d of %d GB entries" % (i + 1, i + len(batch), len(acc_ids)))

            for attempt in range(1, MAX_RETRIES + 1):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="gb", retmode="text")
                    records = handle.read().split("\n//\n")
                    handle.close()

                    for acc_id, record in zip(batch, records):
                        record = record.strip()
                        if not record:
                            continue
                        gbFile = os.path.join(outdir, str(acc_id) + ".gb")
                        with open(gbFile, "w") as outfile:
                            outfile.write(record + "\n//\n")
                        results[acc_id] = gbFile

                    time.sleep(self._sleep_interval)  # stay under NCBI rate limit (3 req/sec without API key)
                    break  # success — move to next batch

                except Exception as e:
                    self.log.warning(
                        "Batch %d-%d failed on attempt %d/%d: %s" % (i + 1, i + len(batch), attempt, MAX_RETRIES,
                                                                     str(e)))
                    if attempt < MAX_RETRIES:
                        self.log.info("Retrying in %d seconds..." % RETRY_DELAY)
                        time.sleep(RETRY_DELAY)
                    else:
                        self.log.error(
                            "Batch %d-%d failed after %d attempts. Falling back to one-by-one fetch." % (i + 1,
                                                                                                         i + len(batch),
                                                                                                         MAX_RETRIES))
                        # Fall back: fetch individually so one bad record doesn't lose the whole batch
                        for acc_id in batch:
                            for solo_attempt in range(1, MAX_RETRIES + 1):
                                try:
                                    gbFile = self.fetch_gb_entry(acc_id, outdir)
                                    results[acc_id] = gbFile
                                    time.sleep(self._sleep_interval)
                                    break
                                except Exception as solo_e:
                                    self.log.warning(
                                        "Solo fetch of %s failed attempt %d/%d: %s" % (acc_id, solo_attempt,
                                                                                       MAX_RETRIES, str(solo_e)))
                                    if solo_attempt < MAX_RETRIES:
                                        time.sleep(RETRY_DELAY)
                                    else:
                                        self.log.error(
                                            "Giving up on accession %s after %d attempts." % (acc_id, MAX_RETRIES))

        return results

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
        fid = fetch_pipe.add_fetch({'retmode':'xml'}, dependency=sid, analyzer=parse_pubmed.PubMedAnalyzer())
        a = cond.run(fetch_pipe)
        result = a.get_result()
        if result.size() >= 1:
            articles = result.pubmed_records
        return articles

    def fetch_xml_entries_batch(self, uids, batch_size=100):
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key

        results = []
        uids = list(uids)
        MAX_RETRIES = 3
        RETRY_DELAY = 5

        for i in range(0, len(uids), batch_size):
            batch = uids[i:i + batch_size]
            self.log.debug("Fetching batch of %d UIDs" % len(batch))

            for attempt in range(1, MAX_RETRIES + 1):
                try:
                    handle = Entrez.efetch(db="nucleotide", id=",".join(map(str, batch)), rettype="gb", retmode="xml")
                    out = handle.read()
                    handle.close()
                    root = ET.fromstring(out)
                    results.extend(root.findall("GBSeq"))
                    time.sleep(self._sleep_interval)
                    break

                except Exception as e:
                    self.log.warning(
                        "XML batch %d-%d failed on attempt %d/%d: %s" % (i + 1, i + len(batch), attempt, MAX_RETRIES,
                                                                         str(e)))
                    if attempt < MAX_RETRIES:
                        self.log.info("Retrying in %d seconds..." % RETRY_DELAY)
                        time.sleep(RETRY_DELAY)
                    else:
                        self.log.error("XML batch %d-%d failed after %d attempts. Skipping." % (i + 1, i + len(batch),
                                                                                                MAX_RETRIES))

        return results

    def internet_on(self):
        try:
            urllib.request.urlopen('https://www.google.com', timeout=5)
            return True
        except urllib.request.URLError:
            try:
                urllib.request.urlopen('https://www.wikipedia.org', timeout=5)
                return True
            except urllib.request.URLError:
                return False


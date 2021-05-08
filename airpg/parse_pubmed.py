import entrezpy.base.analyzer
import entrezpy.base.result
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree

class PubMedRecord():

    def __init__(self):
        self.xml = None

class PubMedResult(entrezpy.base.result.EutilsResult):

    def __init__(self, response, request):
        super().__init__(request.eutil, request.query_id, request.db)
        self.pubmed_records = []

    def size(self):
        return len(self.pubmed_records)

    def isEmpty(self):
        if not self.pubmed_records:
            return True
        return False

    def dump(self):
        return {self:{'dump':{'pubmed_records':[x for x in self.pubmed_records],
                                'query_id': self.query_id, 'db':self.db,
                                'eutil':self.function}}}

    def add_pubmed_record(self, pmc_record):
        self.pubmed_records.append(pmc_record)

class PubMedAnalyzer(entrezpy.base.analyzer.EutilsAnalyzer):

    def __init__(self):
        super().__init__()

    def init_result(self, response, request):
        if self.result is None:
            self.result = PubMedResult(response, request)

    def analyze_error(self, response, request):
        print(json.dumps({__name__:{'Response': {'dump' : request.dump(),
                                               'error' : response.getvalue()}}}))

    def analyze_result(self, response, request):
        self.init_result(response, request)
        articleset = ET.fromstring(response.getvalue())
        for article in list(articleset):
            record = PubMedRecord()
            record.xml = ElementTree(article)
            self.result.add_pubmed_record(record)

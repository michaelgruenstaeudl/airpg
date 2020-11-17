import logging
import re

class ArticleMining:

    def __init__(self, logger):
        self.log = logger or logging.getLogger(__name__ + ".ArticleMining")

    def get_abstract_text(self, article):
        '''
        Parses a pubmed article for its abstract text
        '''
        text = ""
        for elem in article.xml.findall(".//AbstractText"):
            text += elem.text + "\n"
        return text

    def get_article_keywords(self, article):
        '''
        Parses a pubmed article for its keywords
        '''
        keywords = []
        for elem in article.xml.findall(".//Keyword"):
            keywords.append(elem.text)
        #self.log.debug("Found %s keywords." % str(len(keywords)))
        return keywords

    def get_species_from_pubmed_article(self, article, ncbi):
        '''
        Parses a pubmed article for species names
        '''
        species = set()
        abstract_text = self.get_abstract_text(article)
        keywords = self.get_article_keywords(article)
        for keyword in keywords:
            abstract_text += keyword + "\n"
        ncbi_query_results = ncbi.get_name_translator(self.construct_species_query(abstract_text))
        for name, id in ncbi_query_results.items():
            # We're only interested in species-rank taxons
            if ncbi.get_rank(id)[id[0]] == "species":
                species.add(name)

        return species

    def get_genera_from_pubmed_article(self, article, ncbi):
        '''
        Parses a pubmed article for genus names
        '''
        genera = set()
        abstract_text = self.get_abstract_text(article)
        keywords = self.get_article_keywords(article)
        for keyword in keywords:
            abstract_text += keyword + "\n"
        regexp = r'[^a-zA-Z0-9.\-:\s]'
        words = re.sub(regexp, '', abstract_text).split()
        ncbi_results = ncbi.get_name_translator(words)
        for name, id in ncbi_results.items():
            if ncbi.get_rank(id)[id[0]] == "genus":
                genera.add(name)

        return genera


    def construct_species_query(self, text):
        '''
        Constructs a string for every two adjacent words in the text and returns them as a list.
        '''
        species_queries = []
        taxonomic_expressions = ["sp.", "x", "var.", "subsp.", "f."]
        regexp = r'[^a-zA-Z0-9.\-:\s]'
        words = text.split()
        for i in range(0,len(words)-2):
            # Since not all species names are simply two words, we need to account for these cases.
            # Below is an attempt to account for the cases that were detected by a brief look at existing names
            # We're also removing all characters that cannot(should not?) occur in taxonomy (e.g. a species name might be in parentheses)
            if words[i+1] in taxonomic_expressions:
                # "sp." is preceded by one word and is usually followed by one word, but we also found a case with two following words
                if words[i+1] == "sp.":
                    if i < len(words)-3:
                        species_queries.append(re.sub(regexp, '', words[i] + " " + words[i+1] + " " + words[i+2]))
                        if i < len(words)-4:
                            species_queries.append(re.sub(regexp, '', words[i] + " " + words[i+1] + " " + words[i+2] + " " + words[i+3]))
                if words[i+1] == "var." or words[i+1] == "f.":
                    # "var." and "f." are preceded by two words and followed by one word
                    if i < len(words)-3 and i > 0:
                        species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2]))
                if words[i+1] == "subsp.":
                    # "subsp." is preceded by two words and is usually followed by one word, but we also found a case with two following words
                    if i < len(words)-3 and i > 0:
                        species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2]))
                        if i < len(words)-4:
                            species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2] + " " + words[i+3]))
                if words[+1] == "x":
                    # "x" is preceded by two words and followed by two words
                    if i < len(words)-4 and i > 0:
                        species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2] + " " + words[i+3]))
            else:
                species_queries.append(re.sub(regexp, '', words[i] + " " + words[i+1]))
        return species_queries

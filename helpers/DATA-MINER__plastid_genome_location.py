"""
SCRIPT_location_miner.py
------------------------
Mine geographic collection locations for plastid genome (plastome) accessions.

Given a tab-separated input table of plastome accessions (with columns ACCESSION
and TITLE), this script attempts to resolve the country of origin for each entry
via a three-tier fallback pipeline:

  1. Local GenBank records (.tar.gz archives): reads the /country or /lat_lon
     qualifier from the GenBank source feature.
  2. PMC full-text articles (fetched via NCBI Entrez): if no location is found
     in the record, the corresponding publication is retrieved from PubMed Central
     and the article body is searched for keywords such as "collected in" / 
     "collected at", followed by ISO 3166 country-name matching.
  3. Coordinate reverse-geocoding: lat/lon strings are resolved to a country name
     via geopy's Nominatim (OpenStreetMap) geocoder.

The output is a plain-text file with one "count country" line per country,
listing how many accessions originated from each location.

Usage:
    python SCRIPT_location_miner.py -i <input.tsv> -o <output.txt> -e <email>
                                    [-r <records_dir>] [-a <articles_dir>]
                                    [-c] [-v]
"""

import argparse
import io
import os.path
import re
import tarfile
import xml.etree.ElementTree as ET
from collections import Counter

import coloredlogs
import entrezpy.base.analyzer
import entrezpy.base.result
import entrezpy.conduit
import iso3166
import logging
import pandas as pd
from Bio import SeqIO
from geopy.extra.rate_limiter import RateLimiter
from geopy.geocoders import Nominatim


# Keywords searched in PMC article body text to locate collection site sentences
KEYWORDS = ["collected in", "collected at"]
# Regex patterns for detecting degree-format coordinates in free text
# TODO: pattern is currently incomplete (_incomplete placeholder)
PATTERNS = [r"\d{1,3}°\d{1,3}´\s?_incomplete"]


# ---------------------------------------------------------------------------
# PMC / entrezpy helper classes
# ---------------------------------------------------------------------------

class PMCRecord:
    """Holds the raw XML string of a single PubMed Central article."""
    def __init__(self):
        self.xml = None


class PMCResult(entrezpy.base.result.EutilsResult):
    """Accumulates PMCRecord objects returned by an entrezpy fetch pipeline."""
    def __init__(self, response, request):
        super().__init__(request.eutil, request.query_id, request.db)
        self.pmc_records = []

    def size(self):
        return len(self.pmc_records)

    def isEmpty(self):
        return not self.pmc_records

    def dump(self):
        return {self: {'dump': {'pmc_records': list(self.pmc_records),
                                'query_id': self.query_id, 'db': self.db,
                                'eutil': self.function}}}

    def add_pmc_record(self, record):
        self.pmc_records.append(record)


class PMCAnalyzer(entrezpy.base.analyzer.EutilsAnalyzer):
    """entrezpy analyzer that parses PMC XML responses into PMCResult objects."""

    def init_result(self, response, request):
        if self.result is None:
            self.result = PMCResult(response, request)

    def analyze_error(self, response, request):
        import json
        print(json.dumps({__name__: {'Response': {'dump': request.dump(),
                                                   'error': response.getvalue()}}}))

    def analyze_result(self, response, request):
        self.init_result(response, request)
        rec = PMCRecord()
        rec.xml = response.getvalue()
        self.result.add_pmc_record(rec)


# ---------------------------------------------------------------------------
# Country resolution helpers
# ---------------------------------------------------------------------------

def get_country_from_string(text, convert_coordinates=False):
    """
    Attempt to extract a country name from a free-text string.

    First (optionally) scans for degree-format coordinates via PATTERNS and
    resolves them with a geocoder. Then scans each whitespace-separated token
    against the ISO 3166 country list. Returns the first match found, or None.
    """
    if convert_coordinates:
        for pattern in PATTERNS:
            match = re.search(pattern, text)
            if match:
                country = get_country_from_coords(match.group(0))
                if country:
                    return country
    for word in text.split():
        try:
            return iso3166.countries.get(word).name
        except Exception:
            pass
    return None


def get_country_from_record(record, reverse_lookup):
    """
    Extract the country of origin from a Biopython SeqRecord (GenBank format).

    Checks the /country qualifier of the source feature first; falls back to
    reverse-geocoding the /lat_lon qualifier if /country is absent.
    Returns the country name string, or None if not found.
    """
    source = record.features[0]  # first feature should always be of type "source"
    if "country" in source.qualifiers:
        # Some records store "Country: Region" — take only the country part
        return source.qualifiers["country"][0].split(":")[0]
    if "lat_lon" in source.qualifiers:
        return get_country_from_coords(source.qualifiers["lat_lon"], reverse_lookup)
    return None


def get_nominatim_reverse_lookup(user_agent):
    """
    Build and return a rate-limited reverse-geocoding function using
    geopy's Nominatim (OpenStreetMap) geocoder (2-second delay between requests).
    """
    geolocator = Nominatim(user_agent=user_agent, timeout=1)
    return RateLimiter(geolocator.reverse, min_delay_seconds=2)


def get_country_from_coords(coords, reverse_lookup):
    """
    Resolve geographic coordinates to a country name via reverse geocoding.
    Returns the ISO 3166 country name, or None if resolution fails.
    """
    location = reverse_lookup(coords)
    if location:
        try:
            return iso3166.countries.get(location.raw["address"]["country_code"]).name
        except Exception:
            pass
    return None


# ---------------------------------------------------------------------------
# XML helpers
# ---------------------------------------------------------------------------

def get_all_text(xml_elem):
    """Concatenate all text content within an XML element and its descendants."""
    return "".join(elem.text for elem in xml_elem.iter() if elem.text)


def get_section(xml_elem, sec_type):
    """
    Find and return the first <sec> element whose sec-type attribute matches
    sec_type, or None if no such element exists.
    """
    for sec in xml_elem.iter("sec"):
        if sec.attrib.get("sec-type") == sec_type:
            return sec
    return None


def get_title(xml_article):
    """
    Extract the article title string from a parsed PMC XML document.
    Expected element path: article > front > article-meta > title-group > article-title
    TODO: nested tags inside article-title (e.g. <italic>) truncate the returned text.
    """
    return xml_article.find("article/front/article-meta/title-group/article-title").text


# ---------------------------------------------------------------------------
# PMC fetch
# ---------------------------------------------------------------------------

def fetch_pmc_xml(mail, title):
    """
    Search PubMed Central for an article matching the given title and return
    its raw XML string.  The title of the top result is verified against the
    query title; returns None on mismatch or if no results are found.
    """
    cond = entrezpy.conduit.Conduit(mail)
    pipeline = cond.new_pipeline()
    sid = pipeline.add_search({'db': 'pmc', 'term': title, 'rettype': 'count', 'sort': 'relevance'})
    pipeline.add_fetch({'retmax': 1, 'retmode': 'xml'}, dependency=sid, analyzer=PMCAnalyzer())
    result = cond.run(pipeline).get_result()
    if result.size() >= 1:
        pmc_xml = result.pmc_records[0].xml
        if get_title(ET.fromstring(pmc_xml)) == title:
            return pmc_xml
        print("title mismatch")
    return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(args):
    log = logging.getLogger(__name__)
    coloredlogs.install(
        fmt='%(asctime)s [%(levelname)s] %(message)s',
        level='DEBUG' if args.verbose else 'INFO',
        logger=log
    )

    log.debug("Initializing geocoordinate lookup...")
    coordinate_lookup = get_nominatim_reverse_lookup("Plastome_location_miner")
    countries = []           # successfully resolved country names
    nothing_in_records = []  # accessions for which no location was found in local records
    nothing_in_local_articles = []  # accessions for which no location was found in PMC articles

    # STEP 1: Read input table (tab-separated, columns must include ACCESSION and TITLE)
    log.debug("Reading input table...")
    pl_table = pd.read_csv(args.infn, sep='\t', index_col=0, encoding="utf-8")

    # STEP 2: Check local GenBank records (.tar.gz) for /country or /lat_lon qualifiers
    log.info("Checking existing records for location information...")
    if args.recordsdir:
        for acc in pl_table["ACCESSION"]:
            record_path = os.path.join(args.recordsdir, acc + ".tar.gz")
            if os.path.exists(record_path):
                log.debug(f"Checking record '{acc}'...")
                # tar.gz archives contain a single GenBank flat file; wrap in TextIOWrapper
                # because tarfile.extractfile returns a binary BufferedReader
                with tarfile.open(record_path, "r:gz") as tar:
                    rec = SeqIO.read(io.TextIOWrapper(tar.extractfile(tar.getmembers()[0])), "genbank")
                country = get_country_from_record(rec, coordinate_lookup)
                if country:
                    log.debug(f"Found country name '{country}'")
                    countries.append(country)
                else:
                    log.debug("No country information found.")
                    nothing_in_records.append(acc)
            elif args.check_online:
                # Placeholder: online GenBank fetch not yet implemented
                log.debug("fetch and read online record - not implemented yet!")
                nothing_in_records.append(acc)
            else:
                log.debug(f"No local record found for accession '{acc}'.")
                nothing_in_records.append(acc)
    else:
        # No records directory supplied: pass all accessions to the article-mining step
        log.info("Records directory parameter not set. Skipping...")
        nothing_in_records.extend(pl_table["ACCESSION"])

    # STEP 3: For accessions with no record-level location, fetch the PMC full-text
    # article and search the body text for collection-site keywords
    log.info("Checking PMC articles for country information...")
    if args.articledir:
        # TODO: implement reading from a local directory of cached PMC articles
        log.info("Not implemented yet.")
    else:
        for acc in nothing_in_records:
            title = pl_table.loc[pl_table["ACCESSION"] == acc, "TITLE"].values[0]
            log.debug(f"Fetching PMC article titled '{title}'")
            raw_xml = fetch_pmc_xml(args.email, title)
            if raw_xml:
                article = ET.fromstring(raw_xml)
                # Search the entire article body (rather than just Materials & Methods)
                # because section tagging is inconsistent across publishers
                section_text = get_all_text(article.find("article/body"))
                country = None
                for keyword in KEYWORDS:
                    pos = section_text.find(keyword)
                    if pos > -1:
                        log.debug(f"'{keyword}' found! Checking following words for country names...")
                        # Grab the 100 characters immediately after the keyword
                        snippet = section_text[pos + len(keyword):pos + len(keyword) + 100]
                        country = get_country_from_string(snippet, False)
                        if country:
                            break
                if country:
                    log.debug(f"Found country '{country}' mentioned in section.")
                    countries.append(country)
                else:
                    log.debug("No country names or coordinates found in section.")
                    nothing_in_local_articles.append(acc)
            else:
                log.debug(f"Unable to fetch PMC article titled '{title}'")
                nothing_in_local_articles.append(acc)

    # STEP 4: Count occurrences per country and write results
    print(countries)
    log.debug(f"Writing results to output file '{args.outfn}'...")
    with open(args.outfn, "w") as out_fh:
        for country, count in Counter(countries).items():
            out_fh.write(f"{count} {country}\n")
    log.info(f"Could not find country information for {len(nothing_in_local_articles)}/{len(pl_table['ACCESSION'])} accessions.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mine location data from pubmed articles")
    parser.add_argument("-i", "--infn", required=True, help="path to input file")
    parser.add_argument("-o", "--outfn", required=True, help="path to output file")
    parser.add_argument("-r", "--recordsdir", help="path to local GenBank records")
    parser.add_argument("-a", "--articledir", help="path to local PMC articles")
    parser.add_argument("-c", "--check_online", action="store_true", default=False, help="Query GenBank for a record if it does not exist locally.")
    parser.add_argument("-e", "--email", required=True, help="email address for GenBank requests")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="enable verbose logging")
    args = parser.parse_args()
    main(args)

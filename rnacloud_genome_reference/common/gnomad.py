import requests
import json
import logging
from dataclasses import dataclass
from typing import List, Optional

logger = logging.getLogger(__name__)

GNOMAD_REFERENCE_GENOME = 'GRCh38'
GNOMAD_VERSION = 'gnomad_r4'

@dataclass
class GnomadFrequency:
    chrom: str
    pos: int
    alt: str
    lof_filter: Optional[str]
    ac: int
    an: int
    hemizygote_count: int
    homozygote_count: int
    filters: List[str]
    filters_count: int = 0

    def __post_init__(self):
        self.filters_count = len(self.filters)

class GnomadProvider:
    def __init__(self, reference_genome: str = GNOMAD_REFERENCE_GENOME, gnomad_version: str = GNOMAD_VERSION):
        logger.info(f"Initializing GnomadProvider with reference genome: {reference_genome}, gnomAD version: {gnomad_version}")
        self.reference_genome = reference_genome
        self.gnomad_version = gnomad_version
    
    def query_gnomad(self, chrom: str, start: int, stop: int) -> List[GnomadFrequency]:
        """Query gnomAD and return a list of GnomadFrequency objects for the given region."""
        graphql_query = """
            query Region($chrom: String!, $start: Int!, $stop: Int!) {{
            region: region(
                chrom: $chrom
                start: $start
                stop: $stop
                reference_genome: {reference_genome}
            ) {{
                variants(dataset: {gnomad_version}) {{
                chrom
                pos
                alt
                lof_filter
                joint {{
                    ac
                    an
                    hemizygote_count
                    homozygote_count
                    filters
                }}
                }}
            }}
            }}
            """.format(reference_genome=self.reference_genome, gnomad_version=self.gnomad_version)

        url = "https://gnomad.broadinstitute.org/api"
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json",
        }
        variables = {
            "chrom": chrom,
            "start": start,
            "stop": stop
        }
        body = {
            "query": graphql_query,
            "variables": variables
        }
        try:
            logger.debug(f"Query: {graphql_query}")
            logger.debug(f"Variables: {variables}")
            response = requests.post(url, headers=headers, data=json.dumps(body), timeout=600)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"HTTP Request failed: {e}")
            raise ValueError(f"HTTP Request failed for region {chrom}:{start}-{stop}") from e
        data = response.json()
        try:
            variants = data["data"]["region"]["variants"]
        except (KeyError, TypeError):
            print(f"Malformed response or no data for region {chrom}:{start}-{stop}")
            raise ValueError(f"Malformed response or no data for region {chrom}:{start}-{stop}")
        results = []
        for var in variants:
            joint = var.get("joint", {})
            results.append(
                GnomadFrequency(
                    chrom=var.get("chrom"),
                    pos=var.get("pos"),
                    alt=var.get("alt"),
                    lof_filter=var.get("lof_filter"),
                    ac=joint.get("ac", 0),
                    an=joint.get("an", 0),
                    hemizygote_count=joint.get("hemizygote_count", 0),
                    homozygote_count=joint.get("homozygote_count", 0),
                    filters=joint.get("filters", []),
                )
            )
        return results


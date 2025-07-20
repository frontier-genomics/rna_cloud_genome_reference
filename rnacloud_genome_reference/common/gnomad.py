import pandas as pd
import requests
import json
import logging
from dataclasses import dataclass
from typing import Any, List, Optional

logger = logging.getLogger(__name__)

GNOMAD_REFERENCE_GENOME = 'GRCh38'
GNOMAD_VERSION = 'gnomad_r4'

@dataclass
class GnomadFrequency:
    chrom: str
    pos: int
    ref: str
    alt: str
    lof_filter: Optional[str]
    ac: int
    an: int
    hemizygote_count: int
    homozygote_count: int
    filters: List[str]
    filters_count: int = 0
    clinvar_variation_id: Optional[str] = None
    clinical_significance: Optional[str] = None
    review_status: Optional[str] = None

    def __post_init__(self):
        self.filters_count = len(self.filters)

class GnomadProvider:
    def __init__(self, reference_genome: str = GNOMAD_REFERENCE_GENOME, gnomad_version: str = GNOMAD_VERSION):
        logger.info(f"Initializing GnomadProvider with reference genome: {reference_genome}, gnomAD version: {gnomad_version}")
        self.reference_genome = reference_genome
        self.gnomad_version = gnomad_version

    @staticmethod
    def _split_ranges(start: int, stop: int, max_range: int = 50000) -> list[tuple[int, int]]:
        """
        Split a genomic range [start, stop] into subranges no larger than max_range.

        Args:
            start: The starting position (inclusive).
            stop: The ending position (inclusive).
            max_range: Maximum width of each subrange.

        Returns:
            A list of (sub_start, sub_stop) tuples.
        """
        ranges: list[tuple[int, int]] = []
        current_start = start
        while current_start <= stop:
            current_stop = min(current_start + max_range - 1, stop)
            ranges.append((current_start, current_stop))
            current_start = current_stop + 1
        return ranges

    @staticmethod
    def _transform_clinvar_variants(data: list[dict[str, Any]]) -> dict[str, Any]:
        """
        Transforms clinvar_variants from a list to a dict keyed by 'variant_id'.

        Args:
            data (Dict[str, Any]): Input data containing a 'clinvar_variants' list.

        Returns:
            Dict[str, Any]: Output data with 'clinvar_variants' as a dict keyed by variant_id.
        """
        transformed = {}

        for clinvar_entry in data:
            transformed[clinvar_entry['variant_id']] = {
                "clinvar_variation_id": clinvar_entry.get("clinvar_variation_id"),
                "clinical_significance": clinvar_entry.get("clinical_significance"),
                "review_status": clinvar_entry.get("review_status")
            }

        return transformed
    
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
                clinvar_variants {{
                    variant_id
                    clinvar_variation_id
                    clinical_significance
                    review_status
                    major_consequence
                    hgvsc
                }}
                variants(dataset: {gnomad_version}) {{
                    variant_id
                    chrom
                    pos
                    ref
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
            logger.error(f"HTTP Request failed for region {chrom}:{start}-{stop} - {e}")
            return []
        data = response.json()
        try:
            variants = data["data"]["region"]["variants"]
            clinvar_variants = GnomadProvider._transform_clinvar_variants(data["data"]["region"]["clinvar_variants"])
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
                    ref=var.get("ref"),
                    alt=var.get("alt"),
                    lof_filter=var.get("lof_filter"),
                    ac=joint.get("ac", 0),
                    an=joint.get("an", 0),
                    hemizygote_count=joint.get("hemizygote_count", 0),
                    homozygote_count=joint.get("homozygote_count", 0),
                    filters=joint.get("filters", []),
                    clinvar_variation_id=clinvar_variants.get(var.get("variant_id"), {}).get("clinvar_variation_id"),
                    clinical_significance=clinvar_variants.get(var.get("variant_id"), {}).get("clinical_significance"),
                    review_status=clinvar_variants.get(var.get("variant_id"), {}).get("review_status")
                )
            )
        return results

    def fetch_gnomad_stats_for_region(self, chrom: str, start: int, end: int, chunk_size: int = 10000) -> pd.DataFrame | None:
        try:
            total_range = end - start + 1
            if total_range > chunk_size:
                sub_ranges = GnomadProvider._split_ranges(start, end, chunk_size)
                logger.info(
                    f"Requested range {chrom}:{start}-{end} (size={total_range}) "
                    f"exceeds {chunk_size}. Splitting into {len(sub_ranges)} sub-queries."
                )
            else:
                sub_ranges = [(start, end)]

            all_variants: list[GnomadFrequency] = []
            for idx, (sub_start, sub_stop) in enumerate(sub_ranges, start=1):
                logger.info(
                    f"Querying gnomAD ({idx}/{len(sub_ranges)}) for region "
                    f"{chrom}:{sub_start}-{sub_stop}"
                )
                variants = self.query_gnomad(chrom, sub_start, sub_stop)
                if variants:
                    all_variants.extend(variants)
                else:
                    logger.warning(f"No gnomAD data returned for sub-range {sub_start}-{sub_stop}")

            if not all_variants:
                logger.warning(f"No variants found in gnomAD for {chrom}:{start}-{end}")
                return None

            logger.info(f"Total variants found for {chrom}:{start}-{end}: {len(all_variants)}")
            return pd.DataFrame(all_variants)

        except Exception as e:
            logger.error(f"Error querying gnomAD for {chrom}:{start}-{end} - {e}")
            raise ValueError(f"Error querying gnomAD for {chrom}:{start}-{end}") from e
import logging
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from gql import gql, Client
from gql.transport.requests import RequestsHTTPTransport
import pytest
from unittest.mock import Mock, patch

logger = logging.getLogger(__name__)

@dataclass
class GenomicRegion:
    """Represents a genomic region with chromosome, start, and stop positions."""
    chrom: str
    start: int
    stop: int
    id: str | None = None  # Optional ID for the region

    def __post_init__(self):
        if self.start > self.stop:
            raise ValueError(f"Start position {self.start} cannot be greater than stop position {self.stop}")


class GenomicRegionQuerier:
    """Handles querying of multiple genomic regions using GraphQL."""
    
    def __init__(self, endpoint_url: str, reference_genome: str = "GRCh38", dataset: str = "gnomad_r4", 
                 headers: Optional[Dict[str, str]] = None):
        """
        Initialize the querier with GraphQL endpoint and default parameters.
        
        Args:
            endpoint_url: The GraphQL endpoint URL
            reference_genome: Reference genome version (default: GRCh38)
            dataset: Dataset to query (default: gnomad_r4)
            headers: Optional headers for requests
        """
        self.endpoint_url = endpoint_url
        self.reference_genome = reference_genome
        self.dataset = dataset
        self.headers = headers or {}
        
        # Initialize GraphQL client with sync transport
        transport = RequestsHTTPTransport(url=endpoint_url, headers=self.headers)
        self.client = Client(transport=transport, fetch_schema_from_transport=False)
        
        logger.info(f"Initialized GenomicRegionQuerier with endpoint: {endpoint_url}")
    
    def query_multiple_regions(self, regions: List[GenomicRegion]) -> List[Dict[str, Any]]:
        """
        Query multiple genomic regions for variants using GraphQL aliases.
        
        Args:
            regions: List of GenomicRegion objects to query
            
        Returns:
            List of dictionaries containing variant data for each region
            
        Raises:
            ValueError: If regions list is empty
            Exception: If GraphQL query fails
        """
        if not regions:
            raise ValueError("Regions list cannot be empty")
        
        logger.info(f"Querying {len(regions)} genomic regions")
        
        # Build the GraphQL query with aliases
        query_parts = []
        for i, region in enumerate(regions):
            logger.debug(f"Adding region {i+1}: {region.chrom}:{region.start}-{region.stop}")
            
            alias = f"region{i+1}"
            query_part = f"""
                {alias}: region(chrom: "{region.chrom}", start: {region.start}, stop: {region.stop}, reference_genome: {self.reference_genome}) {{
                    variants(dataset: {self.dataset}) {{
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
            """
            query_parts.append(query_part)
        
        # Combine all query parts
        full_query = "query MultipleRegions {\n" + "\n".join(query_parts) + "\n}"
        
        logger.debug(f"Generated GraphQL query. No. of query parts: {len(query_parts)}")
        
        try:
            # Execute the query
            query_doc = gql(full_query)
            
            # Use sync execution with session context
            with self.client as session:
                result = session.execute(query_doc)
            
            logger.info(f"Successfully executed query for {len(regions)} regions")
            
            # Extract and organize results
            results = []
            for i, region in enumerate(regions):
                alias = f"region{i+1}"
                region_data = result.get(alias, {})
                
                # Add region metadata to the result
                region_result = {
                    "region": {
                        "id": region.id,
                        "chrom": region.chrom,
                        "start": region.start,
                        "stop": region.stop
                    },
                    "variants": region_data.get("variants", [])
                }
                results.append(region_result)
                
                variant_count = len(region_result["variants"])
                logger.info(f"Region {region.chrom}:{region.start}-{region.stop} returned {variant_count} variants")
            
            return results
            
        except Exception as e:
            logger.error(f"Failed to execute GraphQL query: {str(e)}")
            raise


def query_genomic_regions(
    regions: List[GenomicRegion], 
    endpoint_url: str,
    reference_genome: str = "GRCh38",
    dataset: str = "gnomad_r4",
    headers: Optional[Dict[str, str]] = None
) -> List[Dict[str, Any]]:
    """
    Convenience function to query multiple genomic regions.
    
    Args:
        regions: List of dictionaries with 'chrom', 'start', 'stop' keys
        endpoint_url: The GraphQL endpoint URL
        reference_genome: Reference genome version (default: GRCh38)
        dataset: Dataset to query (default: gnomad_r4)
        headers: Optional headers for requests
        
    Returns:
        List of dictionaries containing variant data for each region
    """   
    # Create querier and execute query
    querier = GenomicRegionQuerier(endpoint_url, reference_genome, dataset, headers)
    return querier.query_multiple_regions(regions)

if __name__ == "__main__":
    # Example usage
    example_regions = [
        GenomicRegion(chrom="11", start=61398269, stop=61398269, id="example_region_1"),
    ]
    
    endpoint = "https://gnomad.broadinstitute.org/api"
    
    try:
        results = query_genomic_regions(example_regions, endpoint)

        for item in results:
            for variant in item['variants']:
                print(f"{item['region']['id']}\t{variant['chrom']}\t{variant['pos']}\t{variant['alt']}\t{variant['lof_filter']}\t{variant['joint']['ac']}\t{variant['joint']['an']}\t{variant['joint']['hemizygote_count']}\t{variant['joint']['homozygote_count']}")
    except Exception as e:
        logger.error(f"Error querying genomic regions: {str(e)}")
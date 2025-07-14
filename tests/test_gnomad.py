import pytest
from unittest.mock import patch, Mock
from rnacloud_genome_reference.common.gnomad import GnomadProvider, GnomadFrequency  # replace 'your_module' with actual filename

@pytest.mark.integration
def test_query_gnomad_real_call():
    """
    This test will hit the real gnomAD GraphQL API.
    """
    provider = GnomadProvider()
    # Use a region with known data (e.g. chromosome 1, a short region)
    chrom = "1"
    start = 1335374
    stop = 1335376
    
    results = provider.query_gnomad(chrom, start, stop)
    
    assert len(results) == 4, f"Expected 4 results, got {len(results)}"

    results_for_1335374 = list(filter(lambda x: x.pos == 1335374, results))
    assert len(results_for_1335374) == 3, f"Expected 3 result for position 1335374, got {len(results_for_1335374)}"

    results_for_1335376 = list(filter(lambda x: x.pos == 1335376, results))
    assert len(results_for_1335376) == 1, f"Expected 1 result for position 1335374, got {len(results_for_1335376)}"

def test_query_gnomad_mocked():
    """
    Unit test: mock requests.post to test logic in isolation.
    """
    provider = GnomadProvider()

    mock_response = Mock()
    mock_response.json.return_value = {
        "data": {
            "region": {
            "variants": [
                {
                "chrom": "1",
                "pos": 1335374,
                "alt": "G",
                "lof_filter": None,
                "joint": {
                    "ac": 1,
                    "an": 152354,
                    "hemizygote_count": 0,
                    "homozygote_count": 0,
                    "filters": []
                }
                },
                {
                "chrom": "1",
                "pos": 1335374,
                "alt": "C",
                "lof_filter": None,
                "joint": {
                    "ac": 1,
                    "an": 152354,
                    "hemizygote_count": 0,
                    "homozygote_count": 0,
                    "filters": []
                }
                },
                {
                "chrom": "1",
                "pos": 1335374,
                "alt": "T",
                "lof_filter": None,
                "joint": {
                    "ac": 0,
                    "an": 152472,
                    "hemizygote_count": 0,
                    "homozygote_count": 0,
                    "filters": [
                    "AC0"
                    ]
                }
                },
                {
                "chrom": "1",
                "pos": 1335376,
                "alt": "G",
                "lof_filter": None,
                "joint": {
                    "ac": 0,
                    "an": 152474,
                    "hemizygote_count": 0,
                    "homozygote_count": 0,
                    "filters": [
                    "AC0"
                    ]
                }
                }
            ]
            }
        }
        }
    mock_response.raise_for_status = Mock()
    
    def mock_post(*args, **kwargs):
        return mock_response
    
    with patch("requests.post", mock_post):
        results = provider.query_gnomad("1", 1335374, 1335376)

        assert len(results) == 4, f"Expected 4 results, got {len(results)}"

        results_for_1335374 = list(filter(lambda x: x.pos == 1335374, results))
        assert len(results_for_1335374) == 3, f"Expected 3 result for position 1335374, got {len(results_for_1335374)}"

        results_for_1335376 = list(filter(lambda x: x.pos == 1335376, results))
        assert len(results_for_1335376) == 1, f"Expected 1 result for position 1335374, got {len(results_for_1335376)}"
import pytest
from rnacloud_genome_reference.splice_site_population_freq.download_gnomad_freq import split_ranges

@ pytest.mark.parametrize(
    "start, stop, max_range, expected",
    [
        # Single base range
        (1, 1, 50000, [(1, 1)]),
        # Exactly max_range
        (1, 50000, 50000, [(1, 50000)]),
        # One over max_range
        (1, 50001, 50000, [(1, 50000), (50001, 50001)]),
        # Multiple full chunks plus remainder
        (100, 160000, 50000, [(100, 50099), (50100, 100099), (100100, 150099), (150100, 160000)]),
        # start == stop edge case at high value
        (123456, 123456, 50000, [(123456, 123456)]),
        # start greater than stop should produce empty list
        (10, 5, 50000, []),
    ]
)
def test_split_ranges(start: int, stop: int, max_range: int, expected: list[tuple[int, int]]) -> None:
    result = split_ranges(start, stop, max_range)
    assert result == expected, f"Expected {expected} but got {result} for range {start}-{stop} with max_range {max_range}"

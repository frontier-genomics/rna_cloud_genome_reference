import pytest

from rnacloud_genome_reference.genome_build.get_target_contigs import GRC_FIXES_QUERY
from rnacloud_genome_reference.genome_build.generate_mask_bed import range_diff, get_grc_mask_regions, get_cen_par_regions

@pytest.mark.parametrize("start1, end1, start2, end2, expected", [
    # Partial overlap both sides
    (1, 10, 3, 7, [(1, 2), (8, 10)]),
    # Touching at boundary left only
    (1, 5, 5, 8, [(1, 4)]),
    # Touching at boundary right only
    (1, 5, -3, 1, [(2, 5)]),
    # Partial overlap left only
    (1, 10, 7, 15, [(1, 6)]),
    # Partial overlap right only
    (1, 10, -5, 3, [(4, 10)]),
    # Full coverage
    (1, 5, 1, 5, []),
    # No overlap before
    (1, 5, 6, 10, None),
    # No overlap after
    (6, 10, 1, 5, None)
])
def test_range_diff_various(start1, end1, start2, end2, expected):
    assert range_diff(start1, end1, start2, end2) == expected

def test_invalid_first_range_raises():
    with pytest.raises(ValueError):
        range_diff(10, 1, 2, 3)

def test_invalid_second_range_raises():
    with pytest.raises(ValueError):
        range_diff(1, 5, 8, 6)

@pytest.fixture
def QUERY_chr1_MU273333v1_fix():
    return """
        (alt_chr_ucsc == 'chr1_MU273333v1_fix' and clinically_relevant_gene == True)
    """

def test_get_grc_mask_regions_chr1_MU273333v1_fix(QUERY_chr1_MU273333v1_fix):
    result = get_grc_mask_regions(assembly_report='tests/fixtures/GCF_000001405.40_GRCh38.p14_assembly_report.txt',
                                  grc_fixes_assessment='tests/fixtures/grc_fixes_assessment.tsv',
                                  gtf='tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.sorted.gtf.gz',
                                  query=QUERY_chr1_MU273333v1_fix)
    
    results_chr1_MU273333v1_fix = list(filter(lambda r: r.chrom == 'chr1_MU273333v1_fix', result))
    assert len(results_chr1_MU273333v1_fix) == 3, f"Expected 3 regions for chr1_MU273333v1_fix, got {len(results_chr1_MU273333v1_fix)}"

    assert results_chr1_MU273333v1_fix[0].chrom == 'chr1_MU273333v1_fix'
    assert results_chr1_MU273333v1_fix[0].start == 1
    assert results_chr1_MU273333v1_fix[0].end == 1401134
    assert results_chr1_MU273333v1_fix[0].name == 'ATP13A2-SDHB-FIX'

    assert results_chr1_MU273333v1_fix[1].chrom == 'chr1_MU273333v1_fix'
    assert results_chr1_MU273333v1_fix[1].start == 1427112
    assert results_chr1_MU273333v1_fix[1].end == 1433915
    assert results_chr1_MU273333v1_fix[1].name == 'ATP13A2-SDHB-FIX'

    assert results_chr1_MU273333v1_fix[2].chrom == 'chr1_MU273333v1_fix'
    assert results_chr1_MU273333v1_fix[2].start == 1469228
    assert results_chr1_MU273333v1_fix[2].end == 1572686
    assert results_chr1_MU273333v1_fix[2].name == 'ATP13A2-SDHB-FIX'
    

def test_get_cen_par_regions():
    result = get_cen_par_regions(cen_par_regions='tests/fixtures/unmasked_cognates_of_masked_CEN_PAR.txt')

    assert len(result) == 49
    assert result[0].chrom == 'chr5'
    assert result[0].start == 47309185
    assert result[0].end == 49591369
    assert result[0].name == 'GJ212203.1-CEN_PAR'
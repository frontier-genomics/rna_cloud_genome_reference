import pytest

from rnacloud_genome_reference.genome_build.get_target_contigs import get_grc_fixes_contigs, get_ucsc_contigs

@pytest.fixture
def assembly_report() -> str:
    return 'tests/fixtures/GCF_000001405.40_GRCh38.p14_assembly_report.txt'

@pytest.fixture
def grc_fixes_assessment() -> str:
    return 'tests/fixtures/grc_fixes_assessment.tsv'

@pytest.fixture
def GRC_FIXES_QUERY() -> str:
    return '''
        (comparison_status == 'Different - Sequences differ' and clinically_relevant_gene == True)
    '''

@pytest.mark.parametrize("refseq_contigs, ucsc_contigs", [
    (["NT_187388.1"], ["chr22_KI270733v1_random"]),
    (["not existant"], [])
])
def test_get_ucsc_contigs(assembly_report: str, refseq_contigs: list[str], ucsc_contigs: list[str]):
    contigs = get_ucsc_contigs(assembly_report=assembly_report,
                               refseq_contigs=refseq_contigs)

    assert contigs == ucsc_contigs


def test_get_grc_fixes_contigs(grc_fixes_assessment: str, GRC_FIXES_QUERY: str):
    contigs = get_grc_fixes_contigs(grc_fixes_assessment=grc_fixes_assessment,
                                    query=GRC_FIXES_QUERY)

    assert isinstance(contigs, list)
    assert len(contigs) == 22, f"Expected 1 contigs, got {len(contigs)}"

    assert "chr13_ML143365v1_fix" in contigs, f"Expected 'chr13_ML143365v1_fix', but not found"
    assert "chr6_KZ208911v1_fix" in contigs, f"Expected 'chr6_KZ208911v1_fix', but not found"

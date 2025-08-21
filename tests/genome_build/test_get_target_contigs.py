import pytest

from rnacloud_genome_reference.genome_build.get_target_contigs import get_grc_fixes_contigs, get_assembly_report_contigs

@pytest.fixture
def assembly_report() -> str:
    return 'tests/fixtures/GCF_000001405.40_GRCh38.p14_assembly_report.txt'

@pytest.fixture
def grc_fixes_assessment() -> str:
    return 'tests/fixtures/grc_fixes_assessment.tsv'

@pytest.fixture
def ASSEMBLY_REPORT_QUERY() -> str:
    return '''
        `Sequence-Role` == 'assembled-molecule' or \
        `RefSeq-Accn` == 'NT_187388.1' or \
        `RefSeq-Accn` == 'NT_167214.1' or \
        `RefSeq-Accn` == 'NT_187633.1'
    '''

@pytest.fixture
def GRC_FIXES_QUERY() -> str:
    return '''
        (comparison_status == 'Different - Sequences differ' and clinically_relevant_gene == True) or \
        (comparison_status == 'Different - Exon numbering is discordant' and clinically_relevant_gene == True) or \
        (comparison_status == 'Not comparable - Partial transcript annotation in GTF file' and clinically_relevant_gene == True and fix_contig_transcript_partial == False) or \
        (comparison_status == 'Different - No. of exons or introns differ' and clinically_relevant_gene == True)
    '''

def test_get_assembly_report_contigs(assembly_report: str, ASSEMBLY_REPORT_QUERY: str):
    contigs = get_assembly_report_contigs(assembly_report=assembly_report,
                                          query=ASSEMBLY_REPORT_QUERY)

    assert isinstance(contigs, list)
    assert len(contigs) == 28

    assert "chr1" in contigs, f"Expected 'chr1', but not found"
    assert "chr2" in contigs, f"Expected 'chr2', but not found"
    assert "chr3" in contigs, f"Expected 'chr3', but not found"
    assert "chr16_KQ090027v1_alt" not in contigs, f"Expected 'chr16_KQ090027v1_alt' to be excluded, but found"

def test_get_grc_fixes_contigs(grc_fixes_assessment: str, GRC_FIXES_QUERY: str):
    contigs = get_grc_fixes_contigs(grc_fixes_assessment=grc_fixes_assessment,
                                    query=GRC_FIXES_QUERY)

    assert isinstance(contigs, list)
    assert len(contigs) == 47, f"Expected 47 contigs, got {len(contigs)}"

    assert "chr13_ML143365v1_fix" in contigs, f"Expected 'chr13_ML143365v1_fix', but not found"
    assert "chr6_KZ208911v1_fix" in contigs, f"Expected 'chr6_KZ208911v1_fix', but not found"
    assert "chr15_KI270905v1_alt" in contigs, f"Expected 'chr15_KI270905v1_alt', but not found"
    assert "chr19_KI270866v1_alt" in contigs, f"Expected 'chr19_KI270866v1_alt', but not found"

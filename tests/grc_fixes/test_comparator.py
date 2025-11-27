from typing import List

from pytest_bdd import scenario, given, scenarios, then, parsers
from rnacloud_genome_reference.grc_fixes.comparator import FeatureComparisonResult, FeatureSequenceHelper, FeatureComparator
from rnacloud_genome_reference.common.gtf import Exon, Feature, Intron
import pytest

scenarios("features/test_comparator.feature")

def make_intron(sequence, intron_no=1):
    """
    Helper to construct a minimal Intron with a given sequence.
    """
    return Intron(
        chromosome="chr1",
        start=100,
        end=200,
        strand="+",
        sequence=sequence,
        intron_no=intron_no,
    )

class TestFeatureSequenceHelper:
    def test_get_seq_for_region(self):
        input = [
            Exon("NC_000001.11", 65419, 65433, '+', None, 1),
        ]

        expected = [
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1),
        ]

        region_sequence_helper = FeatureSequenceHelper("tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.subset.fna.gz")
        response = region_sequence_helper.get_seq_for_feature(input) # type: ignore

        assert response == expected, f"Expected {expected}, but got {response}"

class TestFeatureComparator:
    @pytest.mark.parametrize("primary_features,fix_features,expected", [
        (
            [Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1)],
            [Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1)],
            0
        ),
        (
            [Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1)],
            [Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAT', 1)],
            1
        ),
        (
            [Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1)],
            [],
            -1
        ),
    ])
    def test_compare_sequences(self, primary_features: list[Feature], fix_features: list[Feature], expected: int):
        response = FeatureComparator.compare_sequences(primary_features, fix_features) # type: ignore

        assert response == expected, f"Expected {expected}, but got {response}"


    def test_compare_features_1(self):
        comparator = FeatureComparator("tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.sorted.gtf.gz", 
                                       "tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.subset.fna.gz")
        
        response = comparator.compare_features(
            primary_chromosome="NC_000021.9",
            primary_start=45405165,
            primary_end=45513720,
            fix_chromosome="NW_025791815.1",
            fix_start=1,
            fix_end=189707,
            entrez_gene_id=80781
        )

        assert isinstance(response, dict), "Response should be a dictionary"
        assert response['primary_contig_transcript'] == 'NM_001379500.1'
        assert response['primary_contig_n_exons'] == 42
        assert response['primary_contig_n_introns'] == 41
        assert response['fix_contig_n_exons'] == 42
        assert response['fix_contig_n_introns'] == 41
        assert response['n_exons_equal'] is True
        assert response['n_introns_equal'] is True
        assert response['sequences_unequal_n_exons'] == 1
        assert response['sequences_unequal_n_introns'] == 0
        assert response['splice_sites_unequal_n'] == 0

    def test_compare_features_2(self):
        comparator = FeatureComparator("tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.sorted.gtf.gz", 
                                       "tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.subset.fna.gz")
        
        response = comparator.compare_features(
            primary_chromosome="NC_000009.12",
            primary_start=134178925,
            primary_end=134206688,
            fix_chromosome="NW_021159999.1",
            fix_start=1,
            fix_end=25408,
            entrez_gene_id=124902298
        )

        assert isinstance(response, dict), "Response should be a dictionary"
        assert response['primary_contig_transcript'] == 'XR_007061835.1'
        assert response['primary_contig_n_exons'] == 3
        assert response['primary_contig_n_introns'] == 2
        assert response['fix_contig_n_exons'] == 0
        assert response['fix_contig_n_introns'] == 0
        assert response['n_exons_equal'] is False
        assert response['n_introns_equal'] is False
        assert response['sequences_unequal_n_exons'] == -1
        assert response['sequences_unequal_n_introns'] == -1
        assert response['splice_sites_unequal_n'] == -1

    
    @pytest.mark.parametrize("primary_seqs, fix_seqs, expected",
        [
            # 1) Both empty → no features to compare
            ([], [], 0),

            # 2) Equal single motif → 0 mismatches
            (["GTAG"], ["GTAG"], 0),

            # 3) 5' motif mismatch only → +1
            (["CTAG"], ["GTAG"], 1),

            # 4) 3' motif mismatch only → +1
            (["ATCC"], ["ATAG"], 1),

            # 5) Both 5' and 3' mismatch → +2
            (["CTCC"], ["GTAG"], 2),

            # 6) Multiple introns: only one of the two differs at either end → +1
            (["GTAG", "CTAG"], ["GTAG", "GTAG"], 1),

            # 8) Length mismatch → returns -1 immediately
            (["GTAG"], ["GTAG", "CTAG"], -1),
        ],
    )
    def test_compare_splice_site_motifs_parametrized(self, primary_seqs: List[str], fix_seqs: List[str], expected: int):
        primary = [
            make_intron(seq, intron_no=i + 1) for i, seq in enumerate(primary_seqs)
        ]
        fix = [
            make_intron(seq, intron_no=i + 1) for i, seq in enumerate(fix_seqs)
        ]

        response = FeatureComparator.compare_splice_site_motifs(primary, fix)
        assert response == expected, f"Expected {expected} but got {response} for primary: {primary_seqs} and fix: {fix_seqs}"

    @pytest.mark.parametrize("primary_seqs, fix_seqs, expected", [
        ([
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1),
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 2)
        ],
        [
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1),
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAT', 2)
        ],
        False),
        ([
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1),
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 2)
        ],
        [
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 2),
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAT', 1)
        ],
        True),
        ([
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1),
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 2)
        ],
        [
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 2)
        ],
        None)
    ])
    def test_flag_discordant_exon_numbering(self, primary_seqs: List[Exon], fix_seqs: List[Exon], expected: bool):
        """
        Test the flag_discordant_exon_numbering method with a specific case.
        """
        response = FeatureComparator.flag_discordant_exon_numbering(primary_seqs, fix_seqs)
        assert response == expected, f"Expected {expected} but got {response} for primary: {primary_seqs} and fix: {fix_seqs}"

   
def _parse_bool_or_none(value: str):
    if value == "True":
        return True
    if value == "False":
        return False
    if value == "None":
        return None
    raise ValueError(f"Cannot parse boolean/None from {value!r}")


def _parse_optional_str(value: str):
    return None if value == "None" else value

@given(
    parsers.parse(
        "a feature comparison {primary_contig_transcript}, "
        "{primary_contig_transcript_is_mane_select}, "
        "{primary_contig_transcript_partial}, "
        "{primary_contig_n_exons}, "
        "{primary_contig_n_introns}, "
        "{fix_contig_transcript}, "
        "{fix_contig_transcript_is_mane_select}, "
        "{fix_contig_transcript_partial}, "
        "{fix_contig_n_exons}, "
        "{fix_contig_n_introns}, "
        "{n_exons_equal}, "
        "{n_introns_equal}, "
        "{sequences_unequal_n_exons}, "
        "{sequences_unequal_n_introns}, "
        "{splice_sites_unequal_n}, "
        "{discordant_exon_numbering}"
    ),
    target_fixture="instance",
)
def feature_comparison_instance(
    primary_contig_transcript,
    primary_contig_transcript_is_mane_select,
    primary_contig_transcript_partial,
    primary_contig_n_exons,
    primary_contig_n_introns,
    fix_contig_transcript,
    fix_contig_transcript_is_mane_select,
    fix_contig_transcript_partial,
    fix_contig_n_exons,
    fix_contig_n_introns,
    n_exons_equal,
    n_introns_equal,
    sequences_unequal_n_exons,
    sequences_unequal_n_introns,
    splice_sites_unequal_n,
    discordant_exon_numbering,
) -> FeatureComparisonResult:
    """
    Build FeatureComparisonResult directly from the Examples table columns.
    All values arrive as strings; we coerce to the proper types here.
    """
    return FeatureComparisonResult(
        primary_contig_transcript=_parse_optional_str(primary_contig_transcript),
        primary_contig_transcript_is_mane_select=_parse_bool_or_none(
            primary_contig_transcript_is_mane_select
        ),
        primary_contig_transcript_partial=_parse_bool_or_none(
            primary_contig_transcript_partial
        ),
        primary_contig_n_exons=int(primary_contig_n_exons),
        primary_contig_n_introns=int(primary_contig_n_introns),
        fix_contig_transcript=_parse_optional_str(fix_contig_transcript),
        fix_contig_transcript_is_mane_select=_parse_bool_or_none(
            fix_contig_transcript_is_mane_select
        ),
        fix_contig_transcript_partial=_parse_bool_or_none(
            fix_contig_transcript_partial
        ),
        fix_contig_n_exons=int(fix_contig_n_exons),
        fix_contig_n_introns=int(fix_contig_n_introns),
        n_exons_equal=_parse_bool_or_none(n_exons_equal),
        n_introns_equal=_parse_bool_or_none(n_introns_equal),
        sequences_unequal_n_exons=int(sequences_unequal_n_exons),
        sequences_unequal_n_introns=int(sequences_unequal_n_introns),
        splice_sites_unequal_n=int(splice_sites_unequal_n),
        discordant_exon_numbering=_parse_bool_or_none(discordant_exon_numbering),
    )


@then(parsers.parse('the comparison status should be "{expected_status}"'))
def comparison_status_should_be(instance: FeatureComparisonResult, expected_status: str):
    assert instance.comparison_status == expected_status
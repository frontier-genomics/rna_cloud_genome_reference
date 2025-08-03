from typing import List
from rnacloud_genome_reference.grc_fixes.comparator import FeatureComparisonResult, FeatureSequenceHelper, FeatureComparator
from rnacloud_genome_reference.common.gtf import Exon, Feature, Intron
import pytest

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

    @pytest.mark.parametrize(
        "instance,expected_status",
        [
            # Identical
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript="tx2",
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=False
                ),
                "Identical",
            ),
            # Different - No. of exons or introns differ
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript="tx2",
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=6,
                    fix_contig_n_introns=4,
                    n_exons_equal=False,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=False
                ),
                "Different - No. of exons or introns differ",
            ),
            # Different - Sequences differ
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript="tx2",
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=1,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=False
                ),
                "Different - Sequences differ",
            ),
            # Different - Splice-site sequences differ
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript="tx2",
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=1,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=False
                ),
                "Different - Splice-site sequences differ",
            ),
            # Different - Exon numbering is discordant
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript="tx2",
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=True,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=False
                ),
                "Different - Exon numbering is discordant",
            ),
            # Not comparable: fix_contig_transcript is None
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript=None,
                    fix_contig_transcript_is_mane_select=False,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=False
                ),
                "Not comparable",
            ),
            (
                FeatureComparisonResult(
                    primary_contig_transcript="NT_187651.1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript=None,
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=True,
                    fix_contig_transcript_partial=None
                ),
                "Not comparable",
            ),
            (
                FeatureComparisonResult(
                    primary_contig_transcript="tx1",
                    primary_contig_transcript_is_mane_select=True,
                    primary_contig_n_exons=5,
                    primary_contig_n_introns=4,
                    fix_contig_transcript='tx1',
                    fix_contig_transcript_is_mane_select=True,
                    fix_contig_n_exons=5,
                    fix_contig_n_introns=4,
                    n_exons_equal=True,
                    n_introns_equal=True,
                    sequences_unequal_n_exons=0,
                    sequences_unequal_n_introns=0,
                    splice_sites_unequal_n=0,
                    discordant_exon_numbering=None,
                    primary_contig_transcript_partial=False,
                    fix_contig_transcript_partial=True
                ),
                "Not comparable - Partial transcript annotation in GTF file",
            )
        ]
    )
    def test_comparison_status(self, instance: FeatureComparisonResult, expected_status: str):
        assert instance.comparison_status == expected_status
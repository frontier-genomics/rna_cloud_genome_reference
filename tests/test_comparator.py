from rnacloud_genome_reference.grc_fixes.comparator import FeatureSequenceHelper, FeatureComparator
from rnacloud_genome_reference.gtf import Exon, Feature, Intron
import pytest

class TestFeatureSequenceHelper:
    def test_get_seq_for_region(self):
        input = [
            Exon("NC_000001.11", 65419, 65433, '+', None, 1),
        ]

        expected = [
            Exon("NC_000001.11", 65419, 65433, '+', 'CCCAGATCTCTTCAG', 1),
        ]

        region_sequence_helper = FeatureSequenceHelper("tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.NC_000001.11.NC_000021.9.NW_025791815.1.NW_025791813.1.NW_025791812.1.fna.gz")
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


    def test_compare_features(self):
        comparator = FeatureComparator("tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.sorted.gtf.gz", 
                                       "tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.NC_000001.11.NC_000021.9.NW_025791815.1.NW_025791813.1.NW_025791812.1.fna.gz")
        
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
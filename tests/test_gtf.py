from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.common.gtf import Exon, Intron
import pytest

class TestGTFHandler:
    @pytest.fixture
    def gtf_hander(self):
        return GTFHandler('tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.sorted.gtf.gz')

    @pytest.mark.parametrize("chromosome, start, end, gene_id, mane, expected", [
        ('NC_000001.11', 65419, 71585, 79501, True, 'NM_001005484.2'),
        ('NC_000001.11', 16538986, 16539663, 124903856, False, 'XM_047436898.1'),
        ('NC_000001.11', 1, 2, 79501, True, None)
    ])
    def test_get_transcript_for_gene(self, gtf_hander: GTFHandler, chromosome: str, start: int, end: int, gene_id: int, mane: bool, expected: str):
        response = gtf_hander.get_transcript_for_gene(chromosome, start, end, gene_id, mane)
        assert response == expected

    @pytest.mark.parametrize("chromosome, start, end, transcript_id, expected", [
        ('NC_000001.11', 65419, 71585, "NM_001005484.2", [
            Exon("NC_000001.11", 65419, 65433, '+', None, 1),
            Exon("NC_000001.11", 65520, 65573, '+', None, 2),
            Exon("NC_000001.11", 69037, 71585, '+', None, 3)
        ]),
        ('NC_000011.10', 5225464, 5226930, "NM_000518.5", [
            Exon("NC_000011.10", 5225464, 5225726, '-', None, 3),
            Exon("NC_000011.10", 5226577, 5226799, '-', None, 2),
            Exon("NC_000011.10", 5226930, 5227071, '-', None, 1)
        ]),
        ('NC_000023.11', 22999960, 23003589, "NM_182699.4", [
            Exon("NC_000023.11", 22999960, 23003589, '+', None, 1)
        ])
    ])
    def test_get_exons_for_transcript(self, gtf_hander: GTFHandler, chromosome: str, start: int, end: int, transcript_id: str, expected: list[Exon]):
        response = gtf_hander.get_exons_by_transcript(chromosome, start, end, transcript_id)
        assert response == expected
        
    @pytest.mark.parametrize("exons, expected", [
        ([
            Exon("NC_000001.11", 65419, 65433, '+', None, 1),
            Exon("NC_000001.11", 65520, 65573, '+', None, 2),
            Exon("NC_000001.11", 69037, 71585, '+', None, 3)
        ], [
            Intron("NC_000001.11", 65434, 65519, '+', None, 1),
            Intron("NC_000001.11", 65574, 69036, '+', None, 2)
        ]),
        ([
            Exon("NC_000011.10", 5225464, 5225726, '-', None, 3),
            Exon("NC_000011.10", 5226577, 5226799, '-', None, 2),
            Exon("NC_000011.10", 5226930, 5227071, '-', None, 1)
        ], [
            Intron("NC_000011.10", 5225727, 5226576, '-', None, 2),
            Intron("NC_000011.10", 5226800, 5226929, '-', None, 1)
        ]),
        ([
            Exon("NC_000023.11", 22999960, 23003589, '+', None, 1)
        ], [])
    ])
    def test_derive_introns_from_exons(self, exons, expected):
        response = GTFHandler.derive_introns_from_exons(exons)
        assert response == expected

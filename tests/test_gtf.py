from rnacloud_genome_reference.gtf import GTFHandler, SpliceJunctionPosition
from rnacloud_genome_reference.gtf import Exon, Intron
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


    @pytest.mark.parametrize("chromosome, start, end, transcript_id, expected", [
        ('NW_012132919.1', 23863, 144565, 'NM_130797.4', True),
        ('NC_000007.14', 153748133, 154894285, 'NM_130797.4', False)
    ])
    def test_is_transcript_partial(self, gtf_hander: GTFHandler, chromosome: str, start: int, end: int, transcript_id: str, expected: bool):
        response = gtf_hander.is_transcript_partial(chromosome, start, end, transcript_id)
        assert response == expected, f"Expected {expected} but got {response} for transcript {transcript_id} in {chromosome}:{start}-{end}"

    @pytest.mark.parametrize("chrom, start, end, entrez_gene_id, expected", [
        ('NC_000001.11', 65419, 71585, 79501, [
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 1, 'Donor', 65434),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 1, 'Donor', 65435),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 2, 'Acceptor', 65518),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 2, 'Acceptor', 65519),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 2, 'Donor', 65574),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 2, 'Donor', 65575),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 3, 'Acceptor', 69035),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', 3, 'Acceptor', 69036)
        ]),
        ('NC_000011.10', 5225464, 5227071, 3043, [
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 3, 'Acceptor', 5225727),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 3, 'Acceptor', 5225728),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 2, 'Donor', 5226575),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 2, 'Donor', 5226576),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 2, 'Acceptor', 5226800),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 2, 'Acceptor', 5226801),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 1, 'Donor', 5226928),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', 1, 'Donor', 5226929)
        ]),
        ('NC_000007.14', 73680918, 73683453, 84277, [])

    ])
    def test_obtain_sj_positions(self, gtf_hander: GTFHandler, chrom: str, start: int, end: int, entrez_gene_id: int, expected: list[SpliceJunctionPosition]):
        response = gtf_hander.obtain_sj_positions(chrom, start, end, entrez_gene_id)

        assert len(response) == len(expected), f"Expected {len(expected)} splice junctions, got {len(response)}"
        
        for item in zip(response, expected):
            assert item[0].chrom == item[1].chrom, f"Chromosome mismatch: {item[0].chrom} != {item[1].chrom}"
            assert item[0].transcript == item[1].transcript, f"Transcript mismatch: {item[0].transcript} != {item[1].transcript}"
            assert item[0].exon_no == item[1].exon_no, f"Exon number mismatch: {item[0].exon_no} != {item[1].exon_no}"
            assert item[0].category == item[1].category, f"Category mismatch: {item[0].category} != {item[1].category}"
            assert item[0].pos == item[1].pos, f"Position mismatch: {item[0].pos} != {item[1].pos}"
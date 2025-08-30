from collections import Counter
from rnacloud_genome_reference.common.gtf import GTFHandler, SpliceJunctionPosition
from rnacloud_genome_reference.common.gtf import Exon, Intron
import pytest

class TestGTFHandler:
    @pytest.fixture
    def gtf_hander(self):
        return GTFHandler('tests/fixtures/GCF_000001405.40_GRCh38.p14_genomic.sorted.gtf.gz')
    
    @pytest.fixture
    def gtf_hander_rna_cloud(self):
        return GTFHandler('tests/fixtures/GCF_000001405.40_GRCh38.p14_rna_cloud_chr1_X.gtf.gz')

    @pytest.mark.parametrize("chromosome, start, end, gene_id, mane, expected", [
        ('NC_000001.11', 65419, 71585, 79501, True, 'NM_001005484.2'),
        ('NC_000001.11', 16538986, 16539663, 124903856, False, 'XM_047436898.1'),
        ('NC_000001.11', 1, 2, 79501, True, None)
    ])
    def test_get_transcript_for_gene(self, gtf_hander: GTFHandler, chromosome: str, start: int, end: int, gene_id: int, mane: bool, expected: str):
        response = gtf_hander.get_transcript_for_gene(chromosome, start, end, gene_id, mane)
        assert response == expected

    def test_get_feature_counts(self, gtf_hander_rna_cloud: GTFHandler):
        response = gtf_hander_rna_cloud.get_feature_counts()
        assert response is not None
        assert isinstance(response, Counter)

        assert len(response) == 6
        
        assert 'CDS' in response.keys()
        assert 'exon' in response.keys()
        assert 'gene' in response.keys()
        assert 'start_codon' in response.keys()
        assert 'stop_codon' in response.keys()
        assert 'transcript' in response.keys()

        assert response['CDS'] == 216746
        assert response['exon'] == 266758
        assert response['gene'] == 5781
        assert response['start_codon'] == 17191
        assert response['stop_codon'] == 17169
        assert response['transcript'] == 23050

    def test_get_gene_biotype_counts(self, gtf_hander_rna_cloud: GTFHandler):
        response = gtf_hander_rna_cloud.get_gene_biotype_counts()
        assert response is not None
        assert isinstance(response, Counter)

        assert len(response) == 14

        assert 'lncRNA' in response.keys()
        assert 'miRNA' in response.keys()
        assert 'transcribed_pseudogene' in response.keys()
        assert 'snoRNA' in response.keys()
        assert 'pseudogene' in response.keys()
        assert 'tRNA' in response.keys()
        assert 'snRNA' in response.keys()
        assert 'rRNA' in response.keys()
        assert 'ncRNA' in response.keys()
        assert 'antisense_RNA' in response.keys()
        assert 'misc_RNA' in response.keys()
        assert 'V_segment' in response.keys()
        assert 'ncRNA_pseudogene' in response.keys()

        assert response['protein_coding'] == 2896
        assert response['lncRNA'] == 1960
        assert response['miRNA'] == 275
        assert response['transcribed_pseudogene'] == 180
        assert response['snoRNA'] == 156
        assert response['pseudogene'] == 143
        assert response['tRNA'] == 89
        assert response['snRNA'] == 41
        assert response['rRNA'] == 17
        assert response['ncRNA'] == 15
        assert response['antisense_RNA'] == 5
        assert response['misc_RNA'] == 2
        assert response['V_segment'] == 1
        assert response['ncRNA_pseudogene'] == 1

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
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 1, 'Donor', 65434, 1),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 1, 'Donor', 65435, 2),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 2, 'Acceptor', 65518, -2),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 2, 'Acceptor', 65519, -1),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 2, 'Donor', 65574, 1),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 2, 'Donor', 65575, 2),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 3, 'Acceptor', 69035, -2),
            SpliceJunctionPosition('NC_000001.11', 'NM_001005484.2', True, 3, 'Acceptor', 69036, -1)
        ]),
        ('NC_000011.10', 5225464, 5227071, 3043, [
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 3, 'Acceptor', 5225727, -1),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 3, 'Acceptor', 5225728, -2),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 2, 'Donor', 5226575, 2),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 2, 'Donor', 5226576, 1),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 2, 'Acceptor', 5226800, -1),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 2, 'Acceptor', 5226801, -2),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 1, 'Donor', 5226928, 2),
            SpliceJunctionPosition('NC_000011.10', 'NM_000518.5', True, 1, 'Donor', 5226929, 1)
        ]),
        ('NC_000007.14', 73680918, 73683453, 84277, []),
        ('NC_000010.11', 103389050, 103396475, 84833, [
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 5, 'Acceptor', 103389167, -1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 5, 'Acceptor', 103389168, -2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 4, 'Donor', 103392189, 2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 4, 'Donor', 103392190, 1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 4, 'Acceptor', 103392284, -1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 4, 'Acceptor', 103392285, -2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 3, 'Donor', 103392369, 2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 3, 'Donor', 103392370, 1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 3, 'Acceptor', 103392467, -1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 3, 'Acceptor', 103392468, -2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 2, 'Donor', 103395744, 2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 2, 'Donor', 103395745, 1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 2, 'Acceptor', 103396033, -1),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 2, 'Acceptor', 103396034, -2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 1, 'Donor', 103396407, 2),
            SpliceJunctionPosition('NC_000010.11', 'NM_001206427.2', True, 1, 'Donor', 103396408, 1)
        ])

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

    @pytest.mark.parametrize("chromosome, entrez_gene_id, start, end, strand", [
        ('NW_012132914.1', 65122, 38599, 43422, '+'),
        ('NC_000001.11', 65122, 12857086, 12861909, '+')
    ])
    def test_get_gene_by_entrez_id(self, gtf_hander: GTFHandler, chromosome: str, entrez_gene_id: int, start: int, end: int, strand: str):
        # Test for a valid Entrez Gene ID
        gene = gtf_hander.get_gene_by_entrez_id(chromosome, entrez_gene_id)
        
        if gene is not None:
            assert gene.chromosome == chromosome, f"Expected chromosome {chromosome}, got {gene.chromosome}"
            assert gene.start == start, f"Expected start {start}, got {gene.start}"
            assert gene.end == end, f"Expected end {end}, got {gene.end}"
            assert gene.strand == strand, f"Expected strand {strand}, got {gene.strand}"


    @pytest.mark.parametrize("chromosome, entrez_gene_id", [
        ('NW_012132914.1', 365122)
    ])
    def test_get_invalid_gene_by_entrez_id(self, gtf_hander: GTFHandler, chromosome: str, entrez_gene_id: int):
        # Test for a valid Entrez Gene ID
        gene = gtf_hander.get_gene_by_entrez_id(chromosome, entrez_gene_id)

        assert gene is None
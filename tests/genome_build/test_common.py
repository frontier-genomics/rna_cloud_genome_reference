import pytest

from rnacloud_genome_reference.genome_build.common import subtract_ranges

def test_ideal_case_simple():
    contig = (1, 10)
    genes = [(3, 4), (8, 9)]
    expected = [(1, 2), (5, 7), (10, 10)]
    assert subtract_ranges(contig, genes) == expected

def test_ideal_case_overlapping_genes():
    contig = (1, 15)
    genes = [(3, 5), (4, 7), (10, 12)]
    expected = [(1, 2), (8, 9), (13, 15)]
    assert subtract_ranges(contig, genes) == expected

def test_ideal_case_full_gene_coverage():
    contig = (1, 5)
    genes = [(1, 5)]
    expected = []
    assert subtract_ranges(contig, genes) == expected

def test_ideal_case_no_genes():
    contig = (1, 5)
    genes = []
    expected = [(1, 5)]
    assert subtract_ranges(contig, genes) == expected

def test_gene_out_of_bounds_start():
    contig = (1, 10)
    genes = [(0, 4)]
    with pytest.raises(ValueError, match="out of contig bounds"):
        subtract_ranges(contig, genes)

def test_gene_out_of_bounds_end():
    contig = (1, 10)
    genes = [(3, 11)]
    with pytest.raises(ValueError, match="out of contig bounds"):
        subtract_ranges(contig, genes)

def test_gene_fully_outside_contig():
    contig = (10, 20)
    genes = [(21, 22)]
    with pytest.raises(ValueError, match="out of contig bounds"):
        subtract_ranges(contig, genes)
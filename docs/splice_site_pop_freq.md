# Splice site population frequency (gnomad v4)

## Output schema
* File is TSV formatted
* Has a header

| Field                     | Description                                        |
|---------------------------|----------------------------------------------------|
| chrom                     | Chromosome                                         |
| chrom_refseq              | Chromosome (Refseq)                                |
| pos                       | Genomic position                                   |
| entrez_gene_id            | Entrez Gene ID                                     |
| gene_name                 | HGNC gene nae                                      |
| transcript                | Transcript                                         |
| transcript_is_mane_select | True if transcript is MANE Select, otherwise False |
| exon_no                   | Exon number                                        |
| dist_from_annot           | Distance from exon                                 |
| category                  | Donor or Acceptor                                  |
| ref                       | Reference motif                                    |
| alt                       | Variant motif                                      |
| lof_filter                | gnomad LoF filter                                  |
| ac                        | Allele count                                       |
| an                        | Allele number                                      |
| af                        | Allele frequency                                   |
| hemizygote_count          | No. of individuals hemizygous for this variant     |
| homozygote_count          | No. of individuals homozygous for this variant     |
| clinvar_variation_id      | ClinVar ID (if available)                          |
| clinical_significance     | Clinical significance (if available)               |
| review_status             | ClinVar review status (if available)               |

## Examples

| **chrom** | **chrom_refseq** | **pos**   | **entrez_gene_id** | **gene_name** | **transcript** | **transcript_is_mane_select** | **exon_no** | **dist_from_annot** | **category** | **ref**     | **alt** | **lof_filter** | **ac** | **an**  | **af**      | **hemizygote_count** | **homozygote_count** | **clinvar_variation_id** | **clinical_significance** | **review_status**                                    |
|-----------|------------------|-----------|--------------------|---------------|----------------|-------------------------------|-------------|---------------------|--------------|-------------|---------|----------------|--------|---------|-------------|----------------------|----------------------|--------------------------|---------------------------|------------------------------------------------------|
| 10        | NC_000010.11     | 103392467 | 84833              | ATP5MK        | NM_001206427.2 | TRUE                          | 3           | -1                  | Acceptor     | CTT         | C       | 5UTR_SPLICE    | 7178   | 1573788 | 0.00456097  | 0                    | 157                  |                          |                           |                                                      |
| 11        | NC_000011.10     | 45935741  | 51317              | PHF21A        | NM_001352027.3 | TRUE                          | 18          | -2                  | Acceptor     | TA          | T       |                | 209315 | 517714  | 0.404306239 | 0                    | 35111                | 403296                   | Benign                    | criteria provided, single submitter                  |
| 11        | NC_000011.10     | 45935741  | 51317              | PHF21A        | NM_001352027.3 | TRUE                          | 18          | -2                  | Acceptor     | TAA         | T       |                | 83030  | 517128  | 0.160559861 | 0                    | 44                   | 403295                   | Benign                    | criteria provided, single submitter                  |
| 12        | NC_000012.12     | 52290076  | 3887               | KRT81         | NM_002281.4    | TRUE                          | 2           | 2                   | Donor        | ACTT        | A       |                | 2391   | 95588   | 0.0250136   | 0                    | 137                  |                          |                           |                                                      |
| 12        | NC_000012.12     | 89472277  | 282809             | POC1B         | NM_172240.3    | TRUE                          | 5           | -2                  | Acceptor     | TAGAAAGAAGA | T       |                | 754313 | 1547452 | 0.487454861 | 0                    | 188722               | 677291                   | Benign                    | criteria provided, multiple submitters, no conflicts |

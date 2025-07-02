
# Gene Primary VS FIX contig comparison schema

## Output schema
* File is TSV formatted
* Has a header

| Field                       | Description                                                                                                                                   |
| --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| chr_refseq                  | Primary contig (Refseq accession) for gene                                                                                                    |
| chr_ucsc                    | Primary contig (UCSC accession) for gene                                                                                                      |
| start                       | Primary contig start position for gene                                                                                                        |
| end                         | Primary contig end position for gene                                                                                                          |
| strand                      | Primary contig strand for gene                                                                                                                |
| gene_name                   | Gene name                                                                                                                                     |
| entrez_gene_id              | NCBI Gene ID / Entrez ID                                                                                                                      |
| issue_id                    | GRC Issue ID(s)                                                                                                                               |
| alt_chr_refseq              | Fix contig (Refseq accession)                                                                                                                 |
| alt_chr_ucsc                | Fix contig (UCSC accession)                                                                                                                   |
| alt_scaf_start              | Fix contig start position                                                                                                                     |
| alt_scaf_stop               | Fix contig end position                                                                                                                       |
| primary_contig_transcript   | Primary contig transcript **<sup>1</sup>**                                                                                                    |
| primary_contig_n_exons      | Primary contig - No. of exons                                                                                                                 |
| primary_contig_n_introns    | Primary contig - No. of introns                                                                                                               |
| fix_contig_transcript       | Same as Primary contig transcript. If Fix contig does not have the same transcript annotated, this is blank and this record is not comparable |
| fix_contig_n_exons          | Fix contig - No. of exons                                                                                                                     |
| fix_contig_n_introns        | Fix contig - No. of introns                                                                                                                   |
| n_exons_equal               | True if Primary contig / Fix contig no. of exons are equal, Otherwise False                                                                   |
| n_introns_equal             | True if Primary contig / Fix contig no. of introns are equal, Otherwise False                                                                 |
| sequences_unequal_n_exons   | No. of exons where sequences differ **<sup>2</sup>**                                                                                          |
| sequences_unequal_n_introns | No. of introns where sequences differ **<sup>2</sup>**                                                                                        |
| comparison_status           | Comparison between gene on primary and FIX contig **<sup>3</sup>**                                                                            |
| clinically_relevant_gene    | If Gene is on the clinically relevant list or not                                                                                             |

**<sup>1</sup>** Primary contig transcript
- Preference is given to MANE
- Otherwise a non MANE transcript is selected (if only one non MANE transcript is present for the gene)

**<sup>2</sup>** Fix contig - No. of exons/introns
- If number of exons/introns differ between Primary/Fix, this is set to -1

**<sup>3</sup>** Comparison between gene on primary and FIX contig
- **Not comparable** - When either the fix contig transcript or primary contig transcript is missing (None)
- **Identical** - When all of the following conditions are met:
    - Number of exons are equal between transcripts
    - Number of introns are equal between transcripts
    - No exon sequences differ (sequences_unequal_n_exons = 0)
    - No intron sequences differ (sequences_unequal_n_introns = 0)
- **Different - No. of exons or introns differ** - When the number of exons or number of introns don't match between the two transcripts
- **Different - Sequences differ** - When the exon and intron counts match, but at least one exon sequence or intron sequence differs between the transcripts
- **Different - Unknown reason** - Fallback status for any other cases that don't match the above conditions

## Examples

| **chr_refseq** | **chr_ucsc** | **start** | **end**  | **strand** | **gene_name** | **entrez_gene_id** | **issue_id**          | **alt_chr_refseq** | **alt_chr_ucsc**    | **alt_scaf_start** | **alt_scaf_stop** | **primary_contig_transcript** | **primary_contig_n_exons** | **primary_contig_n_introns** | **fix_contig_transcript** | **fix_contig_n_exons** | **fix_contig_n_introns** | **n_exons_equal** | **n_introns_equal** | **sequences_unequal_n_exons** | **sequences_unequal_n_introns** | **comparison_status**                      | **clinically_relevant_gene** |
|----------------|--------------|-----------|----------|------------|---------------|--------------------|-----------------------|--------------------|---------------------|--------------------|-------------------|-------------------------------|----------------------------|------------------------------|---------------------------|------------------------|--------------------------|-------------------|---------------------|-------------------------------|---------------------------------|--------------------------------------------|-------------------------------|
| NC_000001.11   | chr1         | 12824605  | 12831410 | -          | PRAMEF11      | 440560             | HG-1342;HG-2282       | NW_012132914.1     | chr1_KQ031383v1_fix | 1                  | 467143            | NM_001146344.3                | 4                          | 3                            | NM_001146344.3            | 4                      | 3                        | TRUE              | TRUE                | 0                             | 0                               | Identical                                  | FALSE                         |
| NC_000001.11   | chr1         | 12847377  | 12848720 | -          | HNRNPCL1      | 343069             | HG-1342;HG-2282       | NW_012132914.1     | chr1_KQ031383v1_fix | 1                  | 467143            | NM_001013631.3                | 2                          | 1                            | NM_001013631.3            | 2                      | 1                        | TRUE              | TRUE                | 0                             | 0                               | Identical                                  | FALSE                         |
| NC_000001.11   | chr1         | 16246840  | 16352480 | -          | FBXO42        | 54455              | HG-1343;HG-173;HG-459 | NW_025791756.1     | chr1_MU273333v1_fix | 1                  | 1572686           | NM_018994.3                   | 10                         | 9                            | NM_018994.3               | 10                     | 9                        | TRUE              | TRUE                | 0                             | 1                               | Different - Sequences differ               | FALSE                         |
| NC_000001.11   | chr1         | 16562423  | 16613564 | -          | NBPF1         | 55672              | HG-1343;HG-173;HG-459 | NW_025791756.1     | chr1_MU273333v1_fix | 1                  | 1572686           | NM_001405678.2                | 28                         | 27                           | NM_001405678.2            | 37                     | 36                       | FALSE             | FALSE               | -1                            | -1                              | Different - No. of exons or introns differ | FALSE                         |
| NC_000001.11   | chr1         | 16672506  | 16673512 | +          | LOC124903857  | 124903857          | HG-1343;HG-173;HG-459 | NW_025791756.1     | chr1_MU273333v1_fix | 1                  | 1572686           | XM_047436909.1                | 1                          | 0                            |                           | 0                      | 0                        | FALSE             | TRUE                | -1                            | 0                               | Not comparable                             | FALSE                         |
| NC_000001.11   | chr1         | 16538986  | 16539663 | +          | LOC124903856  | 124903856          | HG-1343;HG-173;HG-459 | NW_025791756.1     | chr1_MU273333v1_fix | 1                  | 1572686           | XM_047436898.1                | 1                          | 0                            |                           | 0                      | 0                        | FALSE             | TRUE                | -1                            | 0                               | Not comparable                             | FALSE                         |

# RNACloud genome reference

## Workflows

### Gene Primary VS FIX contig comparison

#### Logic

```mermaid
flowchart TD
  G[Download Annotation File]
  H[Download Genome Report]
  I[Download Genome Fasta]
  J[Download GRC Fixes]
  K[Download Clinically Relevant Genes]
    
  G --> N[Extract Protein Coding Genes]
  O[Simplify and Annotate GRC Fixes]
  J --> O
  H --> O

  Q[Determing which protein coding genes overlap FIX contigs]
  O --> Q
  N --> Q
  
  R[Compare Primary and FIX Contigs for each mapping]
  Q --> R
  G --> R
  I --> R

  S[Flag Clinically Relevant Genes]
  R --> S
  K --> S

  S --> T[Final output]

  style G fill:#e1f5fe
  style H fill:#e1f5fe
  style I fill:#e1f5fe
  style J fill:#e1f5fe
  style K fill:#e1f5fe
  style T fill:#f3e5f5,stroke:#000,stroke-width:2px
```
#### Output

https://github.com/frontier-genomics/rnacloud_genome_reference/releases

#### Schema

[Output schema](docs/gene_primary_fix_comparison_summary.md)

## Notes

- NT_187633.1 / chr22_KI270879v1_alt (GSTT1)
# RNACloud Genome Reference - Nextflow Workflow

This directory contains a Nextflow workflow implementation of the RNACloud genome reference processing pipeline.

## Overview

The workflow has been rewritten to use [Nextflow](https://www.nextflow.io/) for better workflow orchestration, reproducibility, and scalability. The workflow maintains the same functionality as the original shell script-based approach but provides:

- **Better reproducibility** through containerization and workflow management
- **Improved scalability** with parallel execution where possible
- **Enhanced portability** across different compute environments
- **Built-in resumability** using Nextflow's `--resume` flag
- **Comprehensive reporting** with timeline, resource usage, and execution reports

## Workflow Components

### Main Processes

1. **Download Processes**: Download required reference files
   - Genome annotation (GTF)
   - Genome sequence (FASTA)
   - Assembly report
   - GRC fixes data
   - Clinically relevant genes

2. **GRC Fixes Assessment**: Assess differences between primary and alternative contigs
   - Extract protein coding genes
   - Process GRC fixes data
   - Compare primary vs alternative sequences
   - Flag clinically relevant genes

3. **Splice Site Population Frequency**: Analyze splice site population frequencies
   - Extract splice junction positions
   - Download gnomAD frequency data
   - Combine splice junctions with population frequency data

## Usage

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (>= 22.04.0)
- [Docker](https://www.docker.com/) (if using Docker profile)
- Java 11 or later (required by Nextflow)

### Quick Start

1. **Build the Docker container** (if needed):
```bash
docker build -t rnacloud_runner .
```

2. **Run the complete workflow**:
```bash
./run_nextflow.sh
```

Or manually:
```bash
nextflow run main.nf --profile docker
```

### Execution Profiles

The workflow supports multiple execution profiles:

- **`standard`** (default): Local execution
- **`docker`**: Docker containerized execution
- **`pbs`**: PBS/Torque cluster execution
- **`slurm`**: SLURM cluster execution

To use a specific profile:
```bash
nextflow run main.nf --profile pbs
```

### Configuration

The workflow uses `nextflow.config` for configuration. Key parameters can be modified:

```bash
# Run with custom output directory
nextflow run main.nf --folders.output_dir /path/to/output

# Run with different gnomAD chunk size
nextflow run main.nf --gnomad.chunk_size 50000
```

### Resume Execution

If the workflow is interrupted, you can resume from where it left off:
```bash
nextflow run main.nf --profile docker --resume
```

## Outputs

The workflow produces the same outputs as the original implementation:

- **GRC Fixes Assessment**: `gene_alt_contigs_mapping_clinically_relevant.tsv`
- **Splice Site Population Frequency**: `*_splice_site_pop_freq.tsv`
- **Workflow Reports**: Timeline, execution report, resource usage trace, and DAG visualization

## Directory Structure

```
.
├── main.nf                    # Main workflow file
├── nextflow.config           # Nextflow configuration
├── run_nextflow.sh          # Convenience script to run workflow
├── workflows/               # Workflow modules
│   ├── download.nf          # Download processes
│   ├── grc_fixes.nf         # GRC fixes assessment
│   └── splice_site.nf       # Splice site frequency analysis
└── rnacloud_genome_reference/  # Original Python modules (unchanged)
```

## Migration from Shell Scripts

The Nextflow workflow replaces the following shell scripts:

- `run_grc_fixes.sh` → Nextflow `ASSESS_GRC_FIXES` process
- `run_splice_site_population_freq.sh` → Nextflow `SPLICE_SITE_POPULATION_FREQ` workflow
- `run_download_gnomad_freq.pbs` → Nextflow `DOWNLOAD_GNOMAD_FREQ` process

All Python modules remain unchanged and are reused by the Nextflow processes.

## Troubleshooting

### Common Issues

1. **Docker permission errors**: Ensure Docker is properly configured for your user
2. **Memory issues**: Adjust memory settings in `nextflow.config` process configuration
3. **Network timeouts**: Increase timeout values for download processes

### Monitoring

Nextflow provides built-in monitoring capabilities:

- **Timeline**: Visual timeline of process execution
- **Report**: Detailed execution statistics
- **Trace**: Process-level resource usage
- **DAG**: Workflow dependency graph

These are automatically generated in the output directory when enabled in the configuration.

## Support

For issues related to:
- **Nextflow workflow**: Check the workflow files in this repository
- **Python modules**: Refer to the original module documentation
- **Nextflow platform**: See [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)
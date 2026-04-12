# Getting started

## Basic usage

The simplest way to get FASTQ files from an SRA accession:

```bash
sracha get SRR000001
```

This will:

1. Resolve the accession via the NCBI SDL API
2. Download the SRA file using parallel chunked HTTP
3. Parse the VDB format natively
4. Output compressed FASTQ files (gzipped by default)

Output files: `SRR000001_1.fastq.gz`, `SRR000001_2.fastq.gz`

## Step by step

If you prefer more control, use the individual subcommands:

```bash
# Download only
sracha fetch SRR000001 -O /data/sra/ --progress

# Convert to FASTQ
sracha fastq /data/sra/SRR000001.sra -O /data/fastq/

# Uncompressed output
sracha fastq SRR000001.sra --no-gzip
```

## SRA-lite

SRA-lite files are smaller (4-10x) because they use simplified quality
scores. To prefer SRA-lite downloads:

```bash
sracha get SRR000001 --format sralite
```

Quality scores will be uniform: Q30 for pass-filter reads, Q3 for rejects.

## Split modes

| Mode | Flag | Output |
|------|------|--------|
| split-3 (default) | `--split split-3` | `_1.fastq.gz`, `_2.fastq.gz`, `_0.fastq.gz` |
| split-files | `--split split-files` | `_1.fastq.gz`, `_2.fastq.gz`, ... |
| split-spot | `--split split-spot` | single file |
| interleaved | `--split interleaved` | single file, R1/R2 alternating |

## Performance tuning

```bash
# More download connections (default: 8)
sracha get SRR000001 --connections 12

# More threads for compression (default: all CPUs)
sracha get SRR000001 --threads 16

# Faster compression (lower ratio)
sracha get SRR000001 --gzip-level 1

# No compression at all
sracha get SRR000001 --no-gzip
```

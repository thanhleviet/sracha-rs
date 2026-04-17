# 🌶️ sracha 🌶️

[![Bioconda](https://anaconda.org/bioconda/sracha/badges/version.svg)](https://anaconda.org/bioconda/sracha) [![Anaconda-Server Badge](https://anaconda.org/bioconda/sracha/badges/downloads.svg)](https://anaconda.org/bioconda/sracha)

Fast SRA downloader and FASTQ converter, written in pure Rust.

## Features

- **Fast** -- 3-10x faster than `fasterq-dump` on typical SRA files
- **One command** -- download, convert to FASTQ, and compress
- **Batch input** -- accessions, BioProjects (PRJNA), studies (SRP), or a file via `--accession-list`
- **gzip or zstd output** -- parallel compression, or plain FASTQ
- **FASTA output** -- `--fasta` drops quality scores
- **SRA and SRA-lite** -- full or simplified quality scores
- **Split modes** -- split-3, split-files, split-spot, interleaved
- **Resumable downloads** -- picks up where it left off
- **Stdout streaming** -- `-Z` pipes FASTQ straight into downstream tools
- **Integrity checks** -- MD5 verification on download and decode
- **Platform support** -- Illumina, BGISEQ/DNBSEQ, Element, Ultima, PacBio, Nanopore (legacy 454 and Ion Torrent are not supported)
- **Single static binary** -- no Python, no C dependencies

## Quick start

```bash
# Download, convert, and compress
sracha get SRR28588231

# Download all runs from a BioProject
sracha get PRJNA675068

# Batch download from an accession list
sracha get --accession-list SRR_Acc_List.txt

# Just download
sracha fetch SRR28588231

# Convert a local .sra file
sracha fastq SRR28588231.sra

# Show accession info
sracha info SRR28588231

# Validate a downloaded file
sracha validate SRR28588231.sra
```

## Benchmarks

Local SRA-to-FASTQ conversion (no network), uncompressed output,
8 CPU cores, measured with [hyperfine](https://github.com/sharkdp/hyperfine).

| File | Size | sracha | fasterq-dump | fastq-dump | Speedup vs fasterq-dump |
|:---|---:|---:|---:|---:|---:|
| SRR28588231 | 23 MiB | 0.21 s | 2.39 s | 2.11 s | **11.2x** |
| SRR2584863 | 288 MiB | 2.01 s | 8.05 s | 13.44 s | **4.0x** |

Compression adds minimal overhead -- sracha produces gzipped FASTQ by default
with parallel block compression, so the integrated pipeline
(`sracha get`) is often faster end-to-end than `fasterq-dump` followed by a
separate gzip step.

<details>
<summary>Full hyperfine output</summary>

**SRR28588231 (23 MiB, 66K spots)**

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `sracha` | 213.2 ± 6.2 | 206.6 | 225.9 | 1.00 |
| `fasterq-dump` | 2385.4 ± 738.8 | 2036.4 | 3706.2 | 11.19 ± 3.48 |
| `fastq-dump` | 2109.8 ± 27.6 | 2082.0 | 2154.8 | 9.89 ± 0.32 |

**SRR2584863 (288 MiB)**

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `sracha` | 2.009 ± 0.045 | 1.981 | 2.062 | 1.00 |
| `fasterq-dump` | 8.052 ± 2.347 | 6.288 | 10.716 | 4.01 ± 1.17 |
| `fastq-dump` | 13.440 ± 1.762 | 12.398 | 15.475 | 6.69 ± 0.89 |

**sracha gzip overhead (SRR28588231)**

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `sracha (no compression)` | 254.9 ± 86.1 | 212.1 | 430.0 | 1.00 |
| `sracha (gzip)` | 342.3 ± 11.1 | 327.3 | 362.1 | 1.34 ± 0.46 |

</details>

Benchmarks run with `sracha` v0.1.11, `sra-tools` v3.2.0, on Linux (8 CPUs).
See `validation/benchmark.sh` to reproduce.

## Installation

Install via [Bioconda](https://bioconda.github.io/):

```bash
pixi add --channel bioconda sracha
```

Or download pre-built binaries from the
[releases page](https://github.com/rnabioco/sracha-rs/releases),
or install from source:

```bash
cargo install --git https://github.com/rnabioco/sracha-rs sracha
```

## Documentation

Full CLI reference and usage guide: <https://rnabioco.github.io/sracha-rs/>

## Acknowledgments

sracha builds on the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra),
maintained by the [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/)
at the National Library of Medicine. The SRA and its
[toolchain](https://github.com/ncbi/sra-tools) are public-domain software
developed by U.S. government employees — our tax dollars at work. Special
thanks to Kenneth Durbrow ([@durbrow](https://github.com/durbrow)) and the
SRA Toolkit team for building and maintaining the infrastructure that makes
projects like this possible.

This project wouldn't exist without NCBI's open infrastructure: the
VDB/KAR format, the SDL locate API, EUtils, and public S3 hosting of
sequencing data. sracha aims to make it easier for the community to
build on that foundation.

## License

MIT

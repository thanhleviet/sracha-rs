# 🌶️ sracha 🌶️

[![Anaconda-Server Badge](https://anaconda.org/bioconda/sracha/badges/version.svg)](https://anaconda.org/bioconda/sracha) [![Anaconda-Server Badge](https://anaconda.org/bioconda/sracha/badges/downloads.svg)](https://anaconda.org/bioconda/sracha)

Fast SRA downloader and FASTQ converter, written in pure Rust.

![sracha demo](docs/images/readme.gif)

## Features

- **Fast** -- 4-11x faster than `fasterq-dump` on typical SRA files
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

### Local decode (SRA file on disk → FASTQ)

Uncompressed output, 8 CPU cores, measured with
[hyperfine](https://github.com/sharkdp/hyperfine).

| File | Size | sracha | fasterq-dump | fastq-dump | Speedup vs fasterq-dump |
|:---|---:|---:|---:|---:|---:|
| SRR28588231 | 23 MiB | 0.16 s | 1.85 s | 2.04 s | **11.6x** |
| SRR2584863 | 288 MiB | 1.29 s | 5.86 s | 12.92 s | **4.5x** |
| ERR1018173 | 1.94 GiB | 7.90 s | 34.47 s | -- | **4.4x** |

Compression adds minimal overhead -- sracha produces gzipped FASTQ by default
with parallel block compression, so the integrated pipeline
(`sracha get`) is often faster end-to-end than `fasterq-dump` followed by a
separate gzip step.

### End-to-end (accession → FASTQ, including download)

Download + decode, 5 runs each from a fresh temp dir. NCBI throughput is
the dominant factor here, so the local-decode speedup is amortized.

| Accession | Size | `sracha get` | `prefetch + fasterq-dump` | `prefetch + fastq-dump` |
|:---|---:|---:|---:|---:|
| SRR28588231 | 23 MiB | 13.54 s | 8.13 s | 8.36 s |
| SRR2584863 | 288 MiB | 19.26 s | 17.05 s | 21.04 s |

On the medium accession `sracha get` is slightly faster than
`prefetch + fastq-dump` and trails `prefetch + fasterq-dump` by ~2 s.
Shared-cluster network variance dominates: `sracha get` std dev was
±1.8 s (small) and ±2.8 s (medium), while prefetch runs stayed within
±0.7 s. Best case (min) `sracha get` matches `fasterq-dump` on the
medium accession. See `validation/bench-results/` for raw hyperfine
output.

<details>
<summary>Full hyperfine output</summary>

**SRR28588231 (23 MiB, 66K spots, Illumina paired)**

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `sracha` | 159.5 ± 16.2 | 139.4 | 198.5 | 1.00 |
| `fasterq-dump` | 1854.7 ± 60.7 | 1789.3 | 1953.5 | 11.63 ± 1.24 |
| `fastq-dump` | 2041.2 ± 37.7 | 2002.4 | 2085.0 | 12.80 ± 1.32 |

**SRR2584863 (288 MiB, Illumina paired)**

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `sracha` | 1.291 ± 0.042 | 1.264 | 1.339 | 1.00 |
| `fasterq-dump` | 5.855 ± 0.097 | 5.746 | 5.933 | 4.53 ± 0.16 |
| `fastq-dump` | 12.915 ± 0.028 | 12.887 | 12.943 | 10.00 ± 0.32 |

**ERR1018173 (1.94 GiB, 15.6M spots, Illumina paired, single run)**

| Command | Time [s] |
|:---|---:|
| `sracha` | 7.90 |
| `fasterq-dump` | 34.47 |

**sracha gzip overhead (SRR28588231)**

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `sracha (no compression)` | 158.1 ± 13.4 | 138.7 | 187.3 | 1.00 |
| `sracha (gzip)` | 318.4 ± 7.4 | 309.5 | 331.4 | 2.01 ± 0.18 |

**End-to-end: SRR28588231 (23 MiB) — accession → FASTQ (5 runs)**

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `sracha get` | 13.542 ± 1.757 | 10.401 | 14.378 | 1.67 ± 0.22 |
| `prefetch + fasterq-dump` | 8.128 ± 0.071 | 8.030 | 8.222 | 1.00 |
| `prefetch + fastq-dump` | 8.363 ± 0.094 | 8.249 | 8.485 | 1.03 ± 0.01 |

**End-to-end: SRR2584863 (288 MiB) — accession → FASTQ (5 runs)**

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `sracha get` | 19.264 ± 2.831 | 16.113 | 21.467 | 1.13 ± 0.17 |
| `prefetch + fasterq-dump` | 17.055 ± 0.663 | 16.541 | 18.125 | 1.00 |
| `prefetch + fastq-dump` | 21.039 ± 1.953 | 19.687 | 24.491 | 1.23 ± 0.12 |

</details>

Benchmarks run with `sracha` v0.3.0, `sra-tools` v3.4.1, on Linux
(8 CPUs). Install the reference toolkit with `pixi run install-sratools`
and reproduce with `validation/benchmark.sh`.

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

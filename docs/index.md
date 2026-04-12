# sracha

Fast SRA downloader and FASTQ converter, written in pure Rust.

## Features

- **Parallel downloads** -- chunked HTTP Range requests with multiple connections
- **Native VDB parsing** -- pure Rust, zero C dependencies
- **Integrated pipeline** -- download, convert, and compress in one command
- **Parallel gzip** -- pigz-style block compression via rayon
- **SRA and SRA-lite** -- full quality or simplified quality scores
- **Split modes** -- split-3, split-files, split-spot, interleaved

## Quick start

```bash
# Download, convert, and compress in one shot
sracha get SRR000001

# Just download
sracha fetch SRR000001

# Convert a local .sra file
sracha fastq SRR000001.sra

# Show accession info
sracha info SRR000001
```

## Installation

### From binary releases

Download pre-built binaries from the
[releases page](https://github.com/rnabioco/sracha-rs/releases).

### From source

```bash
cargo install --git https://github.com/rnabioco/sracha-rs sracha
```

### With pixi

```bash
pixi global install sracha
```

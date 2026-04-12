# CLI reference

## sracha get

Download, convert, and compress SRA data in one shot.

```
sracha get [OPTIONS] <ACCESSION>...
```

### Arguments

| Argument | Description |
|----------|-------------|
| `ACCESSION` | One or more SRA run accessions (SRR/ERR/DRR) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-O, --output-dir` | `.` | Output directory |
| `--format` | `sra` | Preferred format: `sra` or `sralite` |
| `--split` | `split-3` | Split mode: `split-3`, `split-files`, `split-spot`, `interleaved` |
| `--no-gzip` | | Disable gzip (compressed by default) |
| `--gzip-level` | `6` | Compression level (1-9) |
| `-t, --threads` | all CPUs | Thread count |
| `--connections` | `8` | HTTP connections per file |
| `--min-read-len` | | Minimum read length filter |
| `--include-technical` | | Include technical reads |
| `-f, --force` | | Overwrite existing files |
| `-p, --progress` | | Show progress bar |

---

## sracha fetch

Download SRA files without conversion.

```
sracha fetch [OPTIONS] <ACCESSION>...
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-O, --output-dir` | `.` | Output directory |
| `--format` | `sra` | Preferred format: `sra` or `sralite` |
| `--connections` | `8` | HTTP connections per file |
| `--validate` | | Verify MD5 after download |
| `-f, --force` | | Overwrite existing files |
| `-p, --progress` | | Show progress bar |

---

## sracha fastq

Convert local SRA files to FASTQ.

```
sracha fastq [OPTIONS] <INPUT>...
```

### Arguments

| Argument | Description |
|----------|-------------|
| `INPUT` | SRA accession(s) or local `.sra` file path(s) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--split` | `split-3` | Split mode |
| `--no-gzip` | | Disable gzip |
| `--gzip-level` | `6` | Compression level (1-9) |
| `-t, --threads` | all CPUs | Thread count |
| `--min-read-len` | | Minimum read length filter |
| `--include-technical` | | Include technical reads |
| `-Z, --stdout` | | Write to stdout |
| `-O, --output-dir` | `.` | Output directory |
| `-f, --force` | | Overwrite existing files |
| `-p, --progress` | | Show progress bar |

---

## sracha info

Show accession metadata.

```
sracha info <ACCESSION>...
```

Displays file sizes, available formats, download mirrors, and quality
information for each accession.

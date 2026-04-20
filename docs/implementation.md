# Implementation details

## Reference-compressed cSRA

Archives with a physical `CMP_READ` column plus sibling
`PRIMARY_ALIGNMENT` + `REFERENCE` tables decode in pure Rust. Reads
are stored as deltas against a reference genome; sracha reimplements
`NCBI:align:seq_restore_read` and `NCBI:align:align_restore_read` to
splice them back together. Output is byte-identical to `fasterq-dump`
on the fixtures we've verified (e.g. `VDB-3418.sra`).

All FASTQ output flags work on cSRA archives: `--split split-3` /
`--split split-files` / `--split interleaved`, `--fasta`, `--zstd` /
`--no-gzip`, `-Z` stdout streaming, and `-t N` parallel decode all go
through the same FASTQ writer as the plain-SRA path. `sracha fetch`
also pulls the `.sra.vdbcache` sidecar when SDL advertises one so the
AlignmentCursor / ReferenceCursor can read from whichever archive
carries their table.

Two cSRA shapes still error out with "decode with fasterq-dump for
now":

- **Externally-referenced REFERENCE tables** — e.g. SRR341578-class
  archives where `CMP_READ` isn't embedded and REFERENCE bases would
  need to be fetched from refseq.
- **Fixed-length SEQUENCE without physical `READ_LEN`** — the virtual
  cursor that synthesizes `READ_LEN`/`READ_TYPE` from the schema is
  not implemented yet.

Archives that merely *label* themselves with an aligned schema (e.g.
bam-load's `ConvertDatabaseToUnmapped` pathway, which renames
`CMP_READ` → `READ` when no alignments were ingested) keep decoding
through the plain physical-column path.

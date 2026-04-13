//! Pipeline orchestration for the `sracha get` command.
//!
//! Implements a full-download approach:
//!
//! 1. **Phase 1 -- Download**: Use the parallel chunked downloader to fetch the
//!    full SRA file to a temporary location.
//! 2. **Phase 2 -- Parse + Output**: Open the downloaded file as a KAR archive,
//!    create a VdbCursor to read SEQUENCE table columns, decode VDB blobs, and
//!    write FASTQ output.
//! 3. **Phase 3 -- Cleanup + Report**: Delete the temp file and print stats.

use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

use flate2::write::GzEncoder;
use crate::download::{DownloadConfig, download_file};
use crate::error::{Error, Result};
use crate::fastq::{FastqConfig, OutputSlot, SplitMode, SpotRecord, format_spot, output_filename};
use crate::sdl::ResolvedAccession;
use crate::vdb::blob;
use crate::vdb::cursor::VdbCursor;
use crate::vdb::kar::KarArchive;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for the get pipeline.
pub struct PipelineConfig {
    /// Directory for output files.
    pub output_dir: PathBuf,
    /// How to split reads across output files.
    pub split_mode: SplitMode,
    /// Whether to gzip-compress output.
    pub gzip: bool,
    /// Gzip compression level (1-9).
    pub gzip_level: u32,
    /// Number of threads for compression (None = all cores).
    pub threads: Option<usize>,
    /// Number of parallel HTTP connections for downloading.
    pub connections: usize,
    /// Skip technical reads.
    pub skip_technical: bool,
    /// Minimum read length filter.
    pub min_read_len: Option<u32>,
    /// Overwrite existing output files.
    pub force: bool,
    /// Show progress indicators.
    pub progress: bool,
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

/// Statistics from a completed pipeline run.
pub struct PipelineStats {
    /// The accession that was processed.
    pub accession: String,
    /// Number of spots (rows) read from the SRA file.
    pub spots_read: u64,
    /// Number of FASTQ reads written (after filtering).
    pub reads_written: u64,
    /// Total bytes downloaded via HTTP Range requests.
    pub bytes_downloaded: u64,
    /// Total size of the full SRA file on the server.
    pub total_sra_size: u64,
    /// Paths of all output files created.
    pub output_files: Vec<PathBuf>,
}

// ---------------------------------------------------------------------------
// Output writer
// ---------------------------------------------------------------------------

/// An output writer that handles both gzip-compressed and plain FASTQ output.
enum OutputWriter {
    Gz(GzEncoder<std::io::BufWriter<std::fs::File>>),
    Plain(std::io::BufWriter<std::fs::File>),
}

impl OutputWriter {
    fn write_all(&mut self, data: &[u8]) -> std::io::Result<()> {
        match self {
            OutputWriter::Gz(w) => w.write_all(data),
            OutputWriter::Plain(w) => w.write_all(data),
        }
    }

    fn finish(self) -> std::io::Result<()> {
        match self {
            OutputWriter::Gz(w) => {
                w.finish()?;
                Ok(())
            }
            OutputWriter::Plain(mut w) => {
                w.flush()?;
                Ok(())
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Mirror selection
// ---------------------------------------------------------------------------

/// Select the best mirror URL for downloading.
///
/// Prefers cloud mirrors (s3, gs) over NCBI on-premises servers because
/// cloud CDNs are typically much faster for parallel chunked downloads.
fn select_mirror(resolved: &ResolvedAccession) -> Result<String> {
    let mirrors = &resolved.sra_file.mirrors;
    if mirrors.is_empty() {
        return Err(Error::Download {
            accession: resolved.accession.clone(),
            message: "no download mirrors available".into(),
        });
    }

    // Prefer cloud mirrors — much faster for parallel downloads.
    // Priority: s3 > gs > sra-ncbi > ncbi > any
    let priority = |s: &str| -> u8 {
        match s {
            "s3" => 0,
            "gs" => 1,
            s if s.contains("sra-ncbi") => 2,
            "ncbi" => 3,
            _ => 4,
        }
    };

    let best = mirrors
        .iter()
        .min_by_key(|m| priority(m.service.as_str()))
        .unwrap();

    tracing::info!(
        "selected mirror: [{}] {}",
        best.service,
        &best.url[..best.url.len().min(80)],
    );

    Ok(best.url.clone())
}

// ---------------------------------------------------------------------------
// Phase 2: Parse VDB + output FASTQ
// ---------------------------------------------------------------------------

/// Decode VDB columns from a local SRA file, format FASTQ, and write to
/// output files.
///
/// This opens the SRA file as a KAR archive, creates a VdbCursor for the
/// SEQUENCE table, bulk-decompresses each column, and iterates through spots
/// to produce FASTQ output.
fn decode_and_write(
    sra_path: &std::path::Path,
    accession: &str,
    config: &PipelineConfig,
    is_lite: bool,
) -> Result<(u64, u64, Vec<PathBuf>)> {
    let file = std::fs::File::open(sra_path)?;
    let mut archive = KarArchive::open(std::io::BufReader::new(file))?;
    let cursor = VdbCursor::open(&mut archive)?;

    // ------------------------------------------------------------------
    // Streaming blob-by-blob decode and FASTQ output.
    //
    // Instead of loading all blobs into memory at once, process one blob
    // index at a time: decode READ, QUALITY, and READ_LEN for that blob,
    // split into reads, write FASTQ, then free the data before moving on.
    //
    // Physical encodings:
    //   READ:     INSDC:2na:packed (2 bits/base, 4 bases/byte)
    //   QUALITY:  zip_encoding of phred scores (raw-deflate-compressed u8 array)
    //   READ_LEN: izip_encoding of u32 values
    //   READ_TYPE: zip_encoding of u8 values
    //   NAME:     zip_encoding of ASCII bytes
    // ------------------------------------------------------------------

    /// Decode a raw blob, stripping envelope/headers/page_map.
    fn decode_raw(raw: &[u8], checksum_type: u8, row_count: u64) -> Result<blob::DecodedBlob> {
        let cs_size: usize = match checksum_type {
            1 => 4,  // CRC32
            2 => 16, // MD5
            _ => 0,
        };
        let effective = if raw.len() > cs_size {
            &raw[..raw.len() - cs_size]
        } else {
            raw
        };
        blob::decode_blob(effective, 0, row_count, 8)
    }

    /// Decode a zip_encoding data section: the blob header tells us the
    /// version. Version 1 = raw deflate, byte-aligned output. Version 2 =
    /// raw deflate with trailing-bits argument. No headers (v1 blob) = the
    /// data is already the raw-deflate stream.
    ///
    /// Returns decompressed bytes, or falls back to raw bytes on failure.
    fn decode_zip_encoding(decoded: &blob::DecodedBlob) -> Vec<u8> {
        use std::io::Read as _;

        let hdr_version = decoded.headers.first().map(|h| h.version).unwrap_or(0);

        // For zip_encoding, the data section is always raw-deflate
        // compressed (windowBits = -15). The header version only affects
        // whether trailing bits need adjustment (version 2).
        if decoded.data.is_empty() {
            return Vec::new();
        }

        // Try raw deflate decompression.
        let mut dec = flate2::read::DeflateDecoder::new(decoded.data.as_slice());
        let mut out = Vec::new();
        if dec.read_to_end(&mut out).is_ok() && !out.is_empty() {
            // Version 2: trim trailing bits if specified in header args.
            if hdr_version == 2 {
                if let Some(trailing_bits) = decoded.headers.first().and_then(|h| h.args.first()) {
                    let total_bits = out.len() as i64 * 8;
                    let actual_bits = total_bits - (8 - trailing_bits);
                    let actual_bytes = ((actual_bits + 7) / 8) as usize;
                    out.truncate(actual_bytes);
                }
            }
            return out;
        }

        // Fallback: try zlib (with header) in case some columns use it.
        let mut dec2 = flate2::read::ZlibDecoder::new(decoded.data.as_slice());
        let mut out2 = Vec::new();
        if dec2.read_to_end(&mut out2).is_ok() && !out2.is_empty() {
            return out2;
        }

        // Last resort: return raw bytes (might already be uncompressed).
        decoded.data.clone()
    }

    let fastq_config = FastqConfig {
        split_mode: config.split_mode,
        skip_technical: config.skip_technical,
        min_read_len: config.min_read_len,
    };

    // Create output directory.
    std::fs::create_dir_all(&config.output_dir)?;

    // Lazily create output writers as we encounter different output slots.
    let mut writers: HashMap<OutputSlot, OutputWriter> = HashMap::new();
    let mut output_files: Vec<PathBuf> = Vec::new();

    let mut spots_read: u64 = 0;
    let mut reads_written: u64 = 0;

    // Use the READ column as the primary blob source.
    let read_cs = cursor.read_col().meta().checksum_type;
    let num_blobs = cursor.read_col().blob_count();

    tracing::info!(
        "{accession}: streaming decode of {num_blobs} blobs",
    );

    for blob_idx in 0..num_blobs {
        // --------------------------------------------------------------
        // Decode READ blob -> 2na -> ASCII bases.
        // --------------------------------------------------------------
        let read_blob = &cursor.read_col().blobs()[blob_idx];
        let read_raw = cursor.read_col().read_raw_blob_for_row(read_blob.start_id)?;
        let read_decoded = decode_raw(&read_raw, read_cs, read_blob.id_range as u64)?;
        let total_bits = read_decoded.data.len() * 8;
        let adjust = read_decoded.adjust as usize;
        let actual_bases = (total_bits.saturating_sub(adjust)) / 2;
        let read_data = crate::vdb::encoding::unpack_2na(&read_decoded.data, actual_bases);
        let read_page_map = read_decoded.page_map;
        drop(read_raw);

        // --------------------------------------------------------------
        // Decode QUALITY blob (same index).
        // zip_encoding: raw-deflate compressed phred scores.
        // The blob header version tells the unzip variant:
        //   v1 = raw deflate, byte-aligned
        //   v2 = raw deflate with trailing-bits arg
        // After decompression, bytes are numeric phred scores (0-40).
        // --------------------------------------------------------------
        let quality_data: Vec<u8> = if let Some(qcol) = cursor.quality_col() {
            if blob_idx < qcol.blob_count() {
                let qblob = &qcol.blobs()[blob_idx];
                let qcs = qcol.meta().checksum_type;
                let qraw = qcol.read_raw_blob_for_row(qblob.start_id)?;
                let qdecoded = decode_raw(&qraw, qcs, qblob.id_range as u64)?;
                drop(qraw);
                decode_zip_encoding(&qdecoded)
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        // Convert quality to phred+33 ASCII if needed.
        // If quality data is absent or all-zero, treat as SRA-lite.
        let quality_is_empty = quality_data.is_empty()
            || quality_data.iter().take(1000).all(|&b| b == 0);
        let quality_all: Vec<u8> = if !quality_is_empty {
            // Check if values are already ASCII (>= 33) or raw phred (0-40).
            let looks_like_ascii = quality_data.iter().take(100).all(|&b| b >= 33);
            if looks_like_ascii && quality_data.len() == read_data.len() {
                // Already phred+33 ASCII (e.g. SRR000001 454 data).
                quality_data
            } else {
                // Raw phred scores: convert to phred+33 ASCII.
                crate::vdb::encoding::phred_to_ascii(&quality_data)
            }
        } else {
            // SRA-lite: synthetic quality.
            crate::vdb::encoding::sra_lite_quality(read_data.len(), !is_lite)
        };

        // --------------------------------------------------------------
        // Decode READ_LEN blob (same index) -> irzip/izip -> u32 lengths.
        // The page map's row_length tells us reads_per_spot.
        // --------------------------------------------------------------
        let (read_lengths, reads_per_spot): (Vec<u32>, usize) = if let Some(rlcol) = cursor.read_len_col() {
            if blob_idx < rlcol.blob_count() {
                let rlblob = &rlcol.blobs()[blob_idx];
                let rlcs = rlcol.meta().checksum_type;
                let rlraw = rlcol.read_raw_blob_for_row(rlblob.start_id)?;
                let rldecoded = decode_raw(&rlraw, rlcs, rlblob.id_range as u64)?;
                drop(rlraw);

                // The page map row_length tells us how many u32 values
                // per VDB row (i.e. reads per spot).
                let rps = rldecoded.page_map.as_ref()
                    .and_then(|pm| pm.lengths.first().copied())
                    .unwrap_or(1) as usize;

                let hdr_version = rldecoded.headers.first().map(|h| h.version).unwrap_or(0);

                let decoded_ints = if hdr_version >= 1 {
                    let hdr = &rldecoded.headers[0];
                    let planes = hdr.ops.first().copied().unwrap_or(0xFF);
                    let min = hdr.args.first().copied().unwrap_or(0);
                    let slope = hdr.args.get(1).copied().unwrap_or(0);
                    let num_elems = (hdr.osize as u32) / 4;
                    blob::irzip_decode(&rldecoded.data, 32, num_elems, min, slope, planes)
                        .unwrap_or_default()
                } else {
                    let num_elems = rldecoded.row_length
                        .unwrap_or_else(|| (rldecoded.data.len() as u64 * 8) / 32)
                        as u32;
                    blob::izip_decode(&rldecoded.data, 32, num_elems).unwrap_or_default()
                };

                // Expand via page map if present.
                let rl_bytes = if let Some(ref pm) = rldecoded.page_map {
                    let elem_bytes = 4usize; // u32
                    let row_length = pm.lengths.first().copied().unwrap_or(1) as usize;
                    let entry_bytes = row_length * elem_bytes;

                    if !pm.data_runs.is_empty() && pm.data_runs.len() as u64 >= pm.total_rows() {
                        // data_runs contains data_offset indices (random_access variant).
                        let mut expanded = Vec::with_capacity(pm.data_runs.len() * entry_bytes);
                        for &offset in &pm.data_runs {
                            let start = offset as usize * entry_bytes;
                            let end = start + entry_bytes;
                            if end <= decoded_ints.len() {
                                expanded.extend_from_slice(&decoded_ints[start..end]);
                            }
                        }
                        expanded
                    } else if !pm.data_runs.is_empty() {
                        pm.expand_data_runs_bytes(&decoded_ints, elem_bytes)
                    } else {
                        decoded_ints
                    }
                } else {
                    decoded_ints
                };

                let lengths: Vec<u32> = rl_bytes
                    .chunks_exact(4)
                    .map(|chunk| u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                    .collect();

                if blob_idx == 0 {
                    tracing::info!(
                        "READ_LEN: {} values, reads_per_spot={}, first_10={:?}",
                        lengths.len(), rps,
                        &lengths[..lengths.len().min(10)],
                    );
                }

                (lengths, rps)
            } else {
                (vec![read_data.len() as u32], 1)
            }
        } else if let Some(ref pm) = read_page_map {
            // No READ_LEN column: use the READ blob's page map for row lengths.
            let mut row_lengths = Vec::new();
            for (len, run) in pm.lengths.iter().zip(pm.leng_runs.iter()) {
                for _ in 0..*run {
                    if *len > 100_000 {
                        continue;
                    }
                    row_lengths.push(*len);
                }
            }
            (row_lengths, 1)
        } else {
            // No READ_LEN, no page map: treat entire blob as one read.
            (vec![read_data.len() as u32], 1)
        };

        // --------------------------------------------------------------
        // Decode NAME blob (same index) if NAME column exists.
        // zip_encoding: raw-deflate compressed ASCII names.
        // The page map tells us per-row name lengths.
        // --------------------------------------------------------------
        let spot_names: Option<Vec<Vec<u8>>> = if let Some(ncol) = cursor.name_col() {
            if blob_idx < ncol.blob_count() {
                let nblob = &ncol.blobs()[blob_idx];
                let ncs = ncol.meta().checksum_type;
                let nraw = ncol.read_raw_blob_for_row(nblob.start_id)?;
                let ndecoded = decode_raw(&nraw, ncs, nblob.id_range as u64)?;
                drop(nraw);

                let name_bytes = decode_zip_encoding(&ndecoded);

                // Split the concatenated name bytes into per-row names
                // using the page map row lengths.
                let num_spots = if reads_per_spot > 0 {
                    read_lengths.len() / reads_per_spot
                } else {
                    read_lengths.len()
                };

                if let Some(ref pm) = ndecoded.page_map {
                    // Page map tells us per-row name lengths.
                    let mut names = Vec::with_capacity(num_spots);
                    let mut offset = 0usize;
                    for (len, run) in pm.lengths.iter().zip(pm.leng_runs.iter()) {
                        let name_len = *len as usize;
                        for _ in 0..*run {
                            if offset + name_len <= name_bytes.len() {
                                names.push(name_bytes[offset..offset + name_len].to_vec());
                                offset += name_len;
                            }
                        }
                    }
                    if blob_idx == 0 {
                        tracing::info!(
                            "NAME: {} names decoded, first={:?}",
                            names.len(),
                            names.first().map(|n| String::from_utf8_lossy(n).to_string()),
                        );
                    }
                    Some(names)
                } else if let Some(row_len) = ndecoded.row_length {
                    // v1 blob: fixed row_length per row.
                    let rl = row_len as usize;
                    if rl > 0 {
                        let names: Vec<Vec<u8>> = name_bytes
                            .chunks(rl)
                            .map(|c| c.to_vec())
                            .collect();
                        Some(names)
                    } else {
                        None
                    }
                } else {
                    // No page map and no row_length: try NUL or newline splitting.
                    let delimiter = if name_bytes.contains(&0) { 0u8 } else { b'\n' };
                    let names: Vec<Vec<u8>> = name_bytes
                        .split(|&b| b == delimiter)
                        .filter(|s| !s.is_empty())
                        .map(|s| s.to_vec())
                        .collect();
                    if !names.is_empty() { Some(names) } else { None }
                }
            } else {
                None
            }
        } else {
            None
        };

        // TODO: decode READ_TYPE and READ_FILTER from their physical encodings.
        // For now, assume all reads are biological and pass filter.

        // --------------------------------------------------------------
        // Group read_lengths into spots and write FASTQ.
        //
        // reads_per_spot tells us how many consecutive entries in
        // read_lengths belong to one spot. For single-end data this is 1;
        // for paired-end it is typically 2.
        // --------------------------------------------------------------
        let mut seq_offset: usize = 0;
        let mut qual_offset: usize = 0;
        let mut spot_idx_in_blob: usize = 0;

        let rps = reads_per_spot.max(1);
        let mut rl_cursor = 0usize;

        while rl_cursor + rps <= read_lengths.len() {
            let spot_read_lengths = &read_lengths[rl_cursor..rl_cursor + rps];
            let spot_total_bases: usize = spot_read_lengths.iter().map(|&l| l as usize).sum();

            // Extract concatenated sequence for this spot.
            let seq_end = seq_offset + spot_total_bases;
            if seq_end > read_data.len() {
                tracing::warn!(
                    "blob {blob_idx}, spot {spot_idx_in_blob}: sequence overrun at offset \
                     {seq_offset} + {spot_total_bases} > {}; stopping blob",
                    read_data.len(),
                );
                break;
            }
            let sequence = &read_data[seq_offset..seq_end];
            seq_offset = seq_end;

            // Extract concatenated quality for this spot.
            let quality: Vec<u8> = if quality_is_empty {
                crate::vdb::encoding::sra_lite_quality(spot_total_bases, true)
            } else {
                let qual_end = qual_offset + spot_total_bases;
                if qual_end > quality_all.len() {
                    crate::vdb::encoding::sra_lite_quality(spot_total_bases, true)
                } else {
                    let q = quality_all[qual_offset..qual_end].to_vec();
                    qual_offset = qual_end;
                    q
                }
            };

            // Determine spot name: prefer NAME column, fall back to synthetic.
            let name: Vec<u8> = if let Some(ref names) = spot_names {
                if spot_idx_in_blob < names.len() {
                    names[spot_idx_in_blob].clone()
                } else {
                    format!("{}", spots_read as usize + spot_idx_in_blob + 1).into_bytes()
                }
            } else {
                format!("{}", spots_read as usize + spot_idx_in_blob + 1).into_bytes()
            };

            // Build read_types: assume all biological (0) for now.
            let read_types = vec![0u8; rps];
            let read_filter = vec![0u8; rps];

            let spot = SpotRecord {
                name,
                sequence: sequence.to_vec(),
                quality,
                read_lengths: spot_read_lengths.to_vec(),
                read_types,
                read_filter,
                spot_group: Vec::new(),
            };

            let records = format_spot(&spot, accession, &fastq_config);

            for (slot, record) in &records {
                let writer = writers.entry(*slot).or_insert_with(|| {
                    let filename = output_filename(accession, *slot, config.gzip);
                    let path = config.output_dir.join(&filename);
                    output_files.push(path.clone());

                    let file =
                        std::fs::File::create(&path).expect("failed to create output file");
                    let buf = std::io::BufWriter::with_capacity(256 * 1024, file);

                    if config.gzip {
                        OutputWriter::Gz(GzEncoder::new(
                            buf,
                            flate2::Compression::new(config.gzip_level),
                        ))
                    } else {
                        OutputWriter::Plain(buf)
                    }
                });

                writer.write_all(&record.data).map_err(Error::Io)?;
                reads_written += 1;
            }

            rl_cursor += rps;
            spot_idx_in_blob += 1;
        }

        spots_read += spot_idx_in_blob as u64;

        // Log progress every 50 blobs.
        if (blob_idx + 1) % 50 == 0 || blob_idx + 1 == num_blobs {
            tracing::info!(
                "{accession}: decoded {}/{} blobs, {} spots so far",
                blob_idx + 1,
                num_blobs,
                spots_read,
            );
        }

        // Blob data is dropped here, freeing memory before the next blob.
    }

    tracing::info!(
        "{accession}: streaming decode complete -- {} spots, {} reads written",
        spots_read,
        reads_written,
    );

    // Finish all writers.
    for (_, writer) in writers {
        writer.finish().map_err(Error::Io)?;
    }

    Ok((spots_read, reads_written, output_files))
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Statistics from a completed fastq conversion (no download).
pub struct FastqStats {
    /// The accession/label used in FASTQ deflines.
    pub accession: String,
    /// Number of spots (rows) read from the SRA file.
    pub spots_read: u64,
    /// Number of FASTQ reads written (after filtering).
    pub reads_written: u64,
    /// Paths of all output files created.
    pub output_files: Vec<PathBuf>,
}

/// Convert a local SRA file to FASTQ without downloading.
///
/// Opens the SRA file as a KAR archive, creates a VdbCursor, decodes VDB
/// blobs, and writes FASTQ output. The `accession` is used for FASTQ
/// defline naming; if `None`, the filename stem is used.
pub fn run_fastq(
    sra_path: &std::path::Path,
    accession: Option<&str>,
    config: &PipelineConfig,
) -> Result<FastqStats> {
    let acc = accession.map(String::from).unwrap_or_else(|| {
        sra_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string()
    });

    // Detect SRA-lite by checking if the quality column is absent.
    // We pass `false` initially; decode_and_write will handle quality
    // absence gracefully via sra_lite_quality fallback.
    let is_lite = false;

    let (spots_read, reads_written, output_files) =
        decode_and_write(sra_path, &acc, config, is_lite)?;

    Ok(FastqStats {
        accession: acc,
        spots_read,
        reads_written,
        output_files,
    })
}

/// Run the full get pipeline for a single accession.
///
/// This is the full-download approach:
/// 1. Download the complete SRA file to a temporary location
/// 2. Open as KAR archive, parse SEQUENCE table columns
/// 3. Decode VDB blobs, format FASTQ, compress output
/// 4. Delete temporary SRA file
pub async fn run_get(
    resolved: &ResolvedAccession,
    config: &PipelineConfig,
) -> Result<PipelineStats> {
    let accession = &resolved.accession;
    let total_sra_size = resolved.sra_file.size;
    let url = select_mirror(resolved)?;
    let urls = vec![url.clone()];

    tracing::info!("{accession}: starting full download from {url}");

    // --- Phase 1: Download the full SRA file ---
    let temp_filename = format!(".sracha-tmp-{accession}.sra");
    let temp_path = config.output_dir.join(&temp_filename);

    // Ensure output directory exists before downloading.
    tokio::fs::create_dir_all(&config.output_dir).await?;

    let dl_config = DownloadConfig {
        connections: config.connections,
        chunk_size: 0, // adaptive
        force: true,   // always overwrite temp file
        validate: false,
        progress: config.progress,
    };

    tracing::info!(
        "{accession}: downloading {} to {}",
        crate::util::format_size(total_sra_size),
        temp_path.display(),
    );

    let dl_result = download_file(
        &urls,
        total_sra_size,
        resolved.sra_file.md5.as_deref(),
        &temp_path,
        &dl_config,
    )
    .await?;

    let bytes_downloaded = dl_result.size;

    tracing::info!(
        "{accession}: download complete ({})",
        crate::util::format_size(bytes_downloaded),
    );

    // --- Phase 2: Parse VDB + output FASTQ ---
    let is_lite = resolved.sra_file.is_lite;
    let accession_owned = accession.clone();
    let temp_path_for_blocking = temp_path.clone();

    let config_for_blocking = PipelineConfig {
        output_dir: config.output_dir.clone(),
        split_mode: config.split_mode,
        gzip: config.gzip,
        gzip_level: config.gzip_level,
        threads: config.threads,
        connections: config.connections,
        skip_technical: config.skip_technical,
        min_read_len: config.min_read_len,
        force: config.force,
        progress: config.progress,
    };

    let (spots_read, reads_written, output_files) = tokio::task::spawn_blocking(move || {
        decode_and_write(
            &temp_path_for_blocking,
            &accession_owned,
            &config_for_blocking,
            is_lite,
        )
    })
    .await
    .map_err(|e| Error::Vdb(format!("decode task panicked: {e}")))??;

    // --- Phase 3: Cleanup ---
    if let Err(e) = tokio::fs::remove_file(&temp_path).await {
        tracing::warn!(
            "{accession}: failed to remove temp file {}: {e}",
            temp_path.display(),
        );
    }

    tracing::info!(
        "{accession}: done -- {spots_read} spots, {reads_written} reads written, \
         {} downloaded",
        crate::util::format_size(bytes_downloaded),
    );

    Ok(PipelineStats {
        accession: accession.clone(),
        spots_read,
        reads_written,
        bytes_downloaded,
        total_sra_size,
        output_files,
    })
}

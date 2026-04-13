use std::path::{Path, PathBuf};

use md5::{Digest, Md5};
use tokio::io::{AsyncReadExt, AsyncSeekExt, AsyncWriteExt};
use tokio::sync::Semaphore;

use crate::error::{Error, Result};

/// Size thresholds for adaptive chunk sizing.
const SMALL_FILE: u64 = 32 * 1024 * 1024; // 32 MiB
const MEDIUM_FILE: u64 = 256 * 1024 * 1024; // 256 MiB
const LARGE_FILE: u64 = 2 * 1024 * 1024 * 1024; // 2 GiB

/// Maximum retry attempts per chunk.
const MAX_RETRIES: u32 = 3;

/// Configuration for parallel chunked downloads.
pub struct DownloadConfig {
    /// Number of parallel connections (default 8).
    pub connections: usize,
    /// Bytes per chunk (0 = adaptive sizing based on file size).
    pub chunk_size: u64,
    /// Overwrite existing files.
    pub force: bool,
    /// Verify MD5 checksum after download.
    pub validate: bool,
    /// Show progress bar.
    pub progress: bool,
}

impl Default for DownloadConfig {
    fn default() -> Self {
        Self {
            connections: 8,
            chunk_size: 0, // adaptive
            force: false,
            validate: true,
            progress: true,
        }
    }
}

/// Result of a successful download.
pub struct DownloadResult {
    /// Path to the downloaded file.
    pub path: PathBuf,
    /// Size in bytes.
    pub size: u64,
    /// MD5 hex digest (if validation was performed).
    pub md5: Option<String>,
}

/// Information gathered from probing the server.
struct ProbeResult {
    /// Whether the server supports HTTP Range requests.
    supports_range: bool,
    /// Content-Length reported by the server.
    content_length: Option<u64>,
}

/// A byte range for one chunk of the download.
#[derive(Debug, Clone, Copy)]
struct ChunkRange {
    start: u64,
    end: u64, // inclusive
}

impl ChunkRange {
    fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}

/// Choose the chunk size based on file size (adaptive) or use the configured value.
fn effective_chunk_size(file_size: u64, config: &DownloadConfig) -> u64 {
    if config.chunk_size > 0 {
        return config.chunk_size;
    }
    // Adaptive sizing.
    if file_size < MEDIUM_FILE {
        8 * 1024 * 1024 // 8 MiB
    } else if file_size < LARGE_FILE {
        16 * 1024 * 1024 // 16 MiB
    } else {
        32 * 1024 * 1024 // 32 MiB
    }
}

/// Divide the file into byte-range chunks.
fn plan_chunks(file_size: u64, chunk_size: u64) -> Vec<ChunkRange> {
    let mut chunks = Vec::new();
    let mut offset = 0u64;
    while offset < file_size {
        let end = std::cmp::min(offset + chunk_size - 1, file_size - 1);
        chunks.push(ChunkRange { start: offset, end });
        offset = end + 1;
    }
    chunks
}

/// Probe the first URL with a HEAD request to check Range support and size.
async fn probe_url(client: &reqwest::Client, url: &str) -> Result<ProbeResult> {
    let resp = client.head(url).send().await?;

    if !resp.status().is_success() {
        return Err(Error::Download {
            accession: String::new(),
            message: format!("HEAD request failed with HTTP {}", resp.status()),
        });
    }

    let supports_range = resp
        .headers()
        .get(reqwest::header::ACCEPT_RANGES)
        .and_then(|v| v.to_str().ok())
        .is_some_and(|v| v.contains("bytes"));

    // reqwest's content_length() may return 0 for HEAD responses when the
    // body is empty. Read the raw header instead.
    let content_length = resp
        .headers()
        .get(reqwest::header::CONTENT_LENGTH)
        .and_then(|v| v.to_str().ok())
        .and_then(|v| v.parse::<u64>().ok());

    Ok(ProbeResult {
        supports_range,
        content_length,
    })
}

/// Download a single chunk with retries and exponential backoff.
async fn download_chunk(
    client: &reqwest::Client,
    url: &str,
    chunk: ChunkRange,
    output_path: &Path,
    progress: Option<&indicatif::ProgressBar>,
) -> Result<()> {
    let mut last_error = None;

    for attempt in 0..MAX_RETRIES {
        if attempt > 0 {
            // Exponential backoff: 1s, 2s, 4s
            let delay = std::time::Duration::from_secs(1 << attempt);
            tracing::warn!(
                "Retrying chunk {}-{} (attempt {}/{}), backoff {:?}",
                chunk.start,
                chunk.end,
                attempt + 1,
                MAX_RETRIES,
                delay,
            );
            tokio::time::sleep(delay).await;
        }

        match try_download_chunk(client, url, chunk, output_path, progress).await {
            Ok(()) => return Ok(()),
            Err(e) => {
                tracing::warn!(
                    "Chunk {}-{} attempt {} failed: {e}",
                    chunk.start,
                    chunk.end,
                    attempt + 1,
                );
                last_error = Some(e);
            }
        }
    }

    Err(last_error.unwrap_or_else(|| Error::Download {
        accession: String::new(),
        message: "chunk download failed with no error captured".into(),
    }))
}

/// Single attempt to download a chunk.
async fn try_download_chunk(
    client: &reqwest::Client,
    url: &str,
    chunk: ChunkRange,
    output_path: &Path,
    progress: Option<&indicatif::ProgressBar>,
) -> Result<()> {
    let range_header = format!("bytes={}-{}", chunk.start, chunk.end);
    let resp = client
        .get(url)
        .header(reqwest::header::RANGE, &range_header)
        .send()
        .await?;

    let status = resp.status();
    if !status.is_success() && status != reqwest::StatusCode::PARTIAL_CONTENT {
        return Err(Error::Download {
            accession: String::new(),
            message: format!("Range request failed with HTTP {status}"),
        });
    }

    // Open the file for writing at the correct offset (pwrite-style).
    let mut file = tokio::fs::OpenOptions::new()
        .write(true)
        .open(output_path)
        .await?;
    file.seek(std::io::SeekFrom::Start(chunk.start)).await?;

    // Stream the response body in pieces.
    let mut bytes_written = 0u64;
    // Use reqwest's chunk() API which doesn't require futures::StreamExt.
    let mut response = resp;
    while let Some(piece) = response.chunk().await? {
        file.write_all(&piece).await?;
        let n = piece.len() as u64;
        bytes_written += n;
        if let Some(pb) = progress {
            pb.inc(n);
        }
    }

    if bytes_written != chunk.len() {
        return Err(Error::Download {
            accession: String::new(),
            message: format!(
                "chunk {}-{}: expected {} bytes, got {bytes_written}",
                chunk.start,
                chunk.end,
                chunk.len(),
            ),
        });
    }

    Ok(())
}

/// Download a file as a single stream (fallback when server lacks Range support).
async fn download_single_stream(
    client: &reqwest::Client,
    url: &str,
    output_path: &Path,
    expected_size: u64,
    progress: Option<&indicatif::ProgressBar>,
) -> Result<()> {
    let mut response = client.get(url).send().await?;

    if !response.status().is_success() {
        return Err(Error::Download {
            accession: String::new(),
            message: format!("GET request failed with HTTP {}", response.status()),
        });
    }

    let mut file = tokio::fs::File::create(output_path).await?;
    let mut bytes_written = 0u64;

    while let Some(piece) = response.chunk().await? {
        file.write_all(&piece).await?;
        let n = piece.len() as u64;
        bytes_written += n;
        if let Some(pb) = progress {
            pb.inc(n);
        }
    }
    file.flush().await?;

    if expected_size > 0 && bytes_written != expected_size {
        return Err(Error::Download {
            accession: String::new(),
            message: format!("size mismatch: expected {expected_size} bytes, got {bytes_written}"),
        });
    }

    Ok(())
}

/// Compute the MD5 hex digest of a file.
async fn compute_md5(path: &Path) -> Result<String> {
    let mut file = tokio::fs::File::open(path).await?;
    let mut hasher = Md5::new();
    let mut buf = vec![0u8; 64 * 1024]; // 64 KiB read buffer

    loop {
        let n = file.read(&mut buf).await?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }

    let digest = hasher.finalize();
    Ok(digest.iter().map(|b| format!("{b:02x}")).collect())
}

/// Create a progress bar for the download.
fn make_progress_bar(total_size: u64) -> indicatif::ProgressBar {
    use std::io::IsTerminal;

    let pb = if std::io::stderr().is_terminal() {
        indicatif::ProgressBar::new(total_size)
    } else {
        let target = indicatif::ProgressDrawTarget::term_like_with_hz(
            Box::new(crate::pipeline::LogTarget),
            1,
        );
        indicatif::ProgressBar::with_draw_target(Some(total_size), target)
    };
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template(
                "[{elapsed_precise}] [{bar:40}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta})",
            )
            .expect("valid progress bar template")
            .progress_chars("=>-"),
    );
    pb
}

/// Download a file from one or more mirror URLs using parallel HTTP Range requests.
///
/// The downloader probes the first URL for Range support, plans byte-range chunks,
/// and downloads them in parallel with a semaphore limiting concurrency. Falls back
/// to a single-stream download when the server does not support Range requests.
pub async fn download_file(
    urls: &[String],
    expected_size: u64,
    expected_md5: Option<&str>,
    output_path: &Path,
    config: &DownloadConfig,
) -> Result<DownloadResult> {
    if urls.is_empty() {
        return Err(Error::Download {
            accession: String::new(),
            message: "no download URLs provided".into(),
        });
    }

    // Check for existing file.
    if output_path.exists() && !config.force {
        return Err(Error::Download {
            accession: String::new(),
            message: format!(
                "output file already exists: {} (use force to overwrite)",
                output_path.display()
            ),
        });
    }

    // Ensure parent directory exists.
    if let Some(parent) = output_path.parent() {
        tokio::fs::create_dir_all(parent).await?;
    }

    let client = reqwest::Client::builder()
        .user_agent(format!("sracha/{}", env!("CARGO_PKG_VERSION")))
        .build()?;

    // Probe the first URL.
    let url = &urls[0];
    tracing::debug!("probing URL: {url}");
    let probe = probe_url(&client, url).await?;
    tracing::debug!(
        "probe result: range={}, content_length={:?}",
        probe.supports_range,
        probe.content_length
    );

    let file_size = probe.content_length.unwrap_or(expected_size);
    if file_size == 0 {
        return Err(Error::Download {
            accession: String::new(),
            message: "cannot determine file size from server or caller".into(),
        });
    }

    // Set up progress bar (wrapped in Arc so it can be shared across tasks).
    let pb: Option<std::sync::Arc<indicatif::ProgressBar>> = if config.progress {
        Some(std::sync::Arc::new(make_progress_bar(file_size)))
    } else {
        None
    };

    let use_parallel = probe.supports_range && file_size >= SMALL_FILE;

    if use_parallel {
        // --- Parallel chunked download ---
        let chunk_size = effective_chunk_size(file_size, config);
        let chunks = plan_chunks(file_size, chunk_size);

        tracing::debug!(
            "Downloading {} in {} chunks ({} each) with {} connections",
            crate::util::format_size(file_size),
            chunks.len(),
            crate::util::format_size(chunk_size),
            config.connections,
        );

        // Pre-allocate the output file to the expected size.
        {
            let file = tokio::fs::File::create(output_path).await?;
            file.set_len(file_size).await?;
        }

        let semaphore = std::sync::Arc::new(Semaphore::new(config.connections));
        let client = std::sync::Arc::new(client);
        let url = std::sync::Arc::new(url.clone());
        let path = std::sync::Arc::new(output_path.to_path_buf());

        let mut handles = Vec::with_capacity(chunks.len());

        for chunk in chunks {
            let sem = semaphore.clone();
            let cli = client.clone();
            let u = url.clone();
            let p = path.clone();
            let pb_clone = pb.clone();

            let handle = tokio::spawn(async move {
                let _permit = sem.acquire().await.map_err(|e| Error::Download {
                    accession: String::new(),
                    message: format!("semaphore acquire failed: {e}"),
                })?;

                download_chunk(&cli, &u, chunk, &p, pb_clone.as_deref()).await
            });

            handles.push(handle);
        }

        // Await all tasks; collect errors.
        let mut first_error: Option<Error> = None;
        for handle in handles {
            match handle.await {
                Ok(Ok(())) => {}
                Ok(Err(e)) => {
                    if first_error.is_none() {
                        first_error = Some(e);
                    }
                }
                Err(join_err) => {
                    if first_error.is_none() {
                        first_error = Some(Error::Download {
                            accession: String::new(),
                            message: format!("task panicked: {join_err}"),
                        });
                    }
                }
            }
        }

        if let Some(e) = first_error {
            // Clean up partial file on failure.
            let _ = tokio::fs::remove_file(output_path).await;
            return Err(e);
        }
    } else {
        // --- Single-stream fallback ---
        tracing::info!(
            "Downloading {} in single stream (range support: {}, size < 32 MiB: {})",
            crate::util::format_size(file_size),
            probe.supports_range,
            file_size < SMALL_FILE,
        );

        download_single_stream(&client, url, output_path, file_size, pb.as_deref()).await?;
    }

    if let Some(pb) = &pb {
        pb.finish_with_message("download complete");
    }

    // Verify MD5 if requested.
    let md5_hex = if let (true, Some(expected)) = (config.validate, expected_md5) {
        let computed = compute_md5(output_path).await?;
        if computed != expected {
            // Clean up the bad file.
            let _ = tokio::fs::remove_file(output_path).await;
            return Err(Error::ChecksumMismatch {
                expected: expected.to_string(),
                actual: computed,
            });
        }
        Some(computed)
    } else if config.validate {
        // Validate requested but no expected MD5 — compute anyway for the result.
        Some(compute_md5(output_path).await?)
    } else {
        None
    };

    let metadata = tokio::fs::metadata(output_path).await?;

    Ok(DownloadResult {
        path: output_path.to_path_buf(),
        size: metadata.len(),
        md5: md5_hex,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plan_chunks_exact() {
        let chunks = plan_chunks(32, 8);
        assert_eq!(chunks.len(), 4);
        assert_eq!(chunks[0].start, 0);
        assert_eq!(chunks[0].end, 7);
        assert_eq!(chunks[1].start, 8);
        assert_eq!(chunks[1].end, 15);
        assert_eq!(chunks[3].start, 24);
        assert_eq!(chunks[3].end, 31);
    }

    #[test]
    fn test_plan_chunks_remainder() {
        let chunks = plan_chunks(30, 8);
        assert_eq!(chunks.len(), 4);
        assert_eq!(chunks[3].start, 24);
        assert_eq!(chunks[3].end, 29);
        assert_eq!(chunks[3].len(), 6);
    }

    #[test]
    fn test_plan_chunks_single() {
        let chunks = plan_chunks(5, 1024);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].start, 0);
        assert_eq!(chunks[0].end, 4);
    }

    #[test]
    fn test_effective_chunk_size_adaptive() {
        let config = DownloadConfig::default();
        // < 256 MiB => 8 MiB chunks
        assert_eq!(
            effective_chunk_size(100 * 1024 * 1024, &config),
            8 * 1024 * 1024
        );
        // < 2 GiB => 16 MiB chunks
        assert_eq!(
            effective_chunk_size(512 * 1024 * 1024, &config),
            16 * 1024 * 1024
        );
        // >= 2 GiB => 32 MiB chunks
        assert_eq!(
            effective_chunk_size(3 * 1024 * 1024 * 1024, &config),
            32 * 1024 * 1024
        );
    }

    #[test]
    fn test_effective_chunk_size_explicit() {
        let config = DownloadConfig {
            chunk_size: 4 * 1024 * 1024,
            ..Default::default()
        };
        assert_eq!(
            effective_chunk_size(1024 * 1024 * 1024, &config),
            4 * 1024 * 1024
        );
    }

    #[test]
    fn test_default_config() {
        let config = DownloadConfig::default();
        assert_eq!(config.connections, 8);
        assert_eq!(config.chunk_size, 0);
        assert!(!config.force);
        assert!(config.validate);
        assert!(config.progress);
    }

    #[test]
    fn test_chunk_len() {
        let chunk = ChunkRange { start: 0, end: 99 };
        assert_eq!(chunk.len(), 100);

        let chunk = ChunkRange {
            start: 100,
            end: 199,
        };
        assert_eq!(chunk.len(), 100);
    }
}

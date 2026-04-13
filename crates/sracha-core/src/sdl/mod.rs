mod response;

pub use response::*;

use crate::error::{Error, Result};

const SDL_URL: &str = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve";

/// A single download mirror for a resolved file.
#[derive(Debug, Clone)]
pub struct ResolvedMirror {
    /// Direct download URL.
    pub url: String,
    /// Service name (e.g. "ncbi", "s3", "gs").
    pub service: String,
}

/// A resolved SRA file with download URLs and metadata.
#[derive(Debug, Clone)]
pub struct ResolvedFile {
    /// Download mirrors with service labels, filtered to freely accessible locations.
    pub mirrors: Vec<ResolvedMirror>,
    /// File size in bytes.
    pub size: u64,
    /// MD5 hex digest, if provided.
    pub md5: Option<String>,
    /// Whether this is an SRA-lite file (quality scores stripped).
    pub is_lite: bool,
}

/// Resolved download information for a single accession.
#[derive(Debug, Clone)]
pub struct ResolvedAccession {
    /// The accession string.
    pub accession: String,
    /// The primary SRA data file.
    pub sra_file: ResolvedFile,
    /// Optional vdbcache companion file.
    pub vdbcache_file: Option<ResolvedFile>,
}

/// Client for the NCBI SDL (Service Discovery Layer) API.
pub struct SdlClient {
    http: reqwest::Client,
}

impl SdlClient {
    pub fn new() -> Self {
        Self {
            http: reqwest::Client::builder()
                .user_agent(format!("sracha/{}", env!("CARGO_PKG_VERSION")))
                .build()
                .expect("failed to build HTTP client"),
        }
    }

    /// Access the inner reqwest client (for pipeline reuse).
    pub fn http_client(&self) -> &reqwest::Client {
        &self.http
    }

    /// Resolve one or more accessions to download locations via the SDL API.
    ///
    /// Uses GET with query parameters: `?acc=SRR000001&acc=SRR000002`
    pub async fn resolve(&self, accessions: &[String]) -> Result<SdlResponse> {
        let mut url = reqwest::Url::parse(SDL_URL).expect("invalid SDL URL");
        for acc in accessions {
            url.query_pairs_mut().append_pair("acc", acc);
        }

        tracing::debug!("SDL request: {url}");

        let resp = self.http.get(url).send().await.map_err(|e| Error::Sdl {
            message: e.to_string(),
        })?;

        if !resp.status().is_success() {
            return Err(Error::Sdl {
                message: format!("HTTP {}", resp.status()),
            });
        }

        let body = resp.text().await.map_err(|e| Error::Sdl {
            message: e.to_string(),
        })?;
        tracing::debug!("SDL response: {body}");

        let parsed: SdlResponse = serde_json::from_str(&body)?;
        Ok(parsed)
    }

    /// Resolve a single accession and return a clean summary with URLs, size, MD5, etc.
    pub async fn resolve_one(&self, accession: &str) -> Result<ResolvedAccession> {
        let response = self.resolve(&[accession.to_string()]).await?;

        let result = response
            .find_result(accession)
            .ok_or_else(|| Error::NotFound(accession.to_string()))?;

        if let Some(status) = result.status
            && status != 200
        {
            let msg = result.message.as_deref().unwrap_or("unknown error");
            return Err(Error::Sdl {
                message: format!("{accession}: SDL status {status} — {msg}"),
            });
        }

        let sra_sdl = result.find_sra_file().ok_or_else(|| Error::Sdl {
            message: format!("no SRA file found for {accession}"),
        })?;

        let sra_file = resolved_file_from_sdl(sra_sdl)?;

        let vdbcache_file = result
            .find_vdbcache_file()
            .and_then(|f| resolved_file_from_sdl(f).ok());

        Ok(ResolvedAccession {
            accession: accession.to_string(),
            sra_file,
            vdbcache_file,
        })
    }
}

/// Extract a [`ResolvedFile`] from an [`SdlFile`] response entry.
fn resolved_file_from_sdl(sdl_file: &SdlFile) -> Result<ResolvedFile> {
    let mirrors: Vec<ResolvedMirror> = sdl_file
        .locations
        .iter()
        .filter(|loc| !loc.ce_required && !loc.pay_required)
        .map(|loc| ResolvedMirror {
            url: loc.link.clone(),
            service: loc.service.clone().unwrap_or_else(|| "unknown".into()),
        })
        .collect();

    let size = sdl_file.size_bytes().unwrap_or(0);
    let md5 = sdl_file.md5.clone();
    let is_lite = sdl_file.noqual;

    Ok(ResolvedFile {
        mirrors,
        size,
        md5,
        is_lite,
    })
}

impl Default for SdlClient {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_sdl_file(
        file_type: &str,
        noqual: bool,
        size: Option<u64>,
        locations: Vec<SdlLocation>,
    ) -> SdlFile {
        SdlFile {
            file_type: file_type.to_string(),
            name: None,
            accession: None,
            object: None,
            size,
            md5: None,
            format: None,
            modification_date: None,
            noqual,
            locations,
        }
    }

    fn make_location(service: &str, ce: bool, pay: bool) -> SdlLocation {
        SdlLocation {
            link: format!("https://{service}.example.com/file"),
            service: Some(service.to_string()),
            region: None,
            expiration_date: None,
            ce_required: ce,
            pay_required: pay,
        }
    }

    #[test]
    fn resolved_file_filters_paid_mirrors() {
        let sdl = make_sdl_file(
            "sra",
            false,
            Some(1000),
            vec![
                make_location("s3", false, false),
                make_location("gs", true, false),   // CE required
                make_location("ncbi", false, true), // pay required
            ],
        );
        let resolved = resolved_file_from_sdl(&sdl).unwrap();
        assert_eq!(resolved.mirrors.len(), 1);
        assert_eq!(resolved.mirrors[0].service, "s3");
    }

    #[test]
    fn resolved_file_keeps_free_mirrors() {
        let sdl = make_sdl_file(
            "sra",
            false,
            Some(5000),
            vec![
                make_location("s3", false, false),
                make_location("ncbi", false, false),
            ],
        );
        let resolved = resolved_file_from_sdl(&sdl).unwrap();
        assert_eq!(resolved.mirrors.len(), 2);
        assert_eq!(resolved.size, 5000);
    }

    #[test]
    fn resolved_file_lite_flag() {
        let sdl = make_sdl_file("sra", true, Some(1000), vec![]);
        let resolved = resolved_file_from_sdl(&sdl).unwrap();
        assert!(resolved.is_lite);
    }
}

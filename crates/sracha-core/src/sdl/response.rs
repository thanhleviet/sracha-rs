use serde::Deserialize;
use serde_json::Value;

/// Top-level response from the NCBI SDL API.
#[derive(Debug, Deserialize)]
pub struct SdlResponse {
    /// API version string (e.g., "2").
    pub version: Option<String>,
    /// Status code (integer in JSON, e.g., 200 or 500).
    #[serde(default, deserialize_with = "deserialize_status")]
    pub status: Option<i64>,
    pub message: Option<String>,
    #[serde(default, rename = "nextToken")]
    pub next_token: Option<String>,
    #[serde(default, rename = "result")]
    pub results: Vec<SdlResult>,
}

impl SdlResponse {
    /// Find the result entry matching a given accession string.
    pub fn find_result(&self, accession: &str) -> Option<&SdlResult> {
        self.results.iter().find(|r| {
            r.query.as_deref() == Some(accession) || r.bundle.as_deref() == Some(accession)
        })
    }
}

/// Per-accession result from SDL.
#[derive(Debug, Deserialize)]
pub struct SdlResult {
    /// The accession that was queried.
    pub query: Option<String>,
    /// Alternate field name used in some responses.
    pub bundle: Option<String>,
    /// Status code (integer, e.g., 200).
    #[serde(default, deserialize_with = "deserialize_status")]
    pub status: Option<i64>,
    #[serde(default, rename = "msg")]
    pub message: Option<String>,
    #[serde(default)]
    pub files: Vec<SdlFile>,
}

impl SdlResult {
    /// Get the accession string from whichever field is populated.
    pub fn accession(&self) -> &str {
        self.query
            .as_deref()
            .or(self.bundle.as_deref())
            .unwrap_or("unknown")
    }

    /// Find the first file entry with type "sra".
    pub fn find_sra_file(&self) -> Option<&SdlFile> {
        self.files.iter().find(|f| f.file_type == "sra")
    }

    /// Find the first file entry with type "vdbcache".
    pub fn find_vdbcache_file(&self) -> Option<&SdlFile> {
        self.files.iter().find(|f| f.file_type == "vdbcache")
    }

    /// Whether this result indicates success.
    pub fn is_ok(&self) -> bool {
        self.status == Some(200)
    }
}

/// A downloadable file for an accession.
#[derive(Debug, Deserialize)]
pub struct SdlFile {
    /// File type: "sra", "vdbcache", etc.
    #[serde(rename = "type")]
    pub file_type: String,
    /// File name.
    pub name: Option<String>,
    /// Accession (sometimes present).
    pub accession: Option<String>,
    /// Object identifier.
    pub object: Option<String>,
    /// File size in bytes (integer in JSON).
    #[serde(default, deserialize_with = "deserialize_size")]
    pub size: Option<u64>,
    /// MD5 hex digest.
    pub md5: Option<String>,
    /// File format.
    pub format: Option<String>,
    /// Modification date.
    #[serde(rename = "modificationDate")]
    pub modification_date: Option<String>,
    /// Whether quality scores are absent (SRA-lite).
    #[serde(default)]
    pub noqual: bool,
    /// Download locations.
    #[serde(default)]
    pub locations: Vec<SdlLocation>,
}

impl SdlFile {
    /// File size in bytes.
    pub fn size_bytes(&self) -> Option<u64> {
        self.size
    }

    /// Whether this is an SRA (not vdbcache or other) file.
    pub fn is_sra(&self) -> bool {
        self.file_type == "sra"
    }
}

/// A specific download location/mirror for a file.
#[derive(Debug, Deserialize)]
pub struct SdlLocation {
    /// Direct download URL.
    pub link: String,
    /// Service name: "ncbi", "s3", "gs".
    pub service: Option<String>,
    /// Geographic region: "be-md", "us-east-1", etc.
    pub region: Option<String>,
    /// URL expiration date (for cloud presigned URLs).
    #[serde(rename = "expirationDate")]
    pub expiration_date: Option<String>,
    /// Whether Compute Environment token is required.
    #[serde(default, rename = "ceRequired")]
    pub ce_required: bool,
    /// Whether payment is required.
    #[serde(default, rename = "payRequired")]
    pub pay_required: bool,
}

// ---------------------------------------------------------------------------
// Flexible deserializers (SDL returns integers but some fields could be strings)
// ---------------------------------------------------------------------------

/// Deserialize a status field that may be an integer or a string.
fn deserialize_status<'de, D>(deserializer: D) -> std::result::Result<Option<i64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let v = Option::<Value>::deserialize(deserializer)?;
    match v {
        None => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_i64()),
        Some(Value::String(s)) => Ok(s.parse().ok()),
        _ => Ok(None),
    }
}

/// Deserialize a size field that may be an integer or a string.
fn deserialize_size<'de, D>(deserializer: D) -> std::result::Result<Option<u64>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let v = Option::<Value>::deserialize(deserializer)?;
    match v {
        None => Ok(None),
        Some(Value::Number(n)) => Ok(n.as_u64()),
        Some(Value::String(s)) => Ok(s.parse().ok()),
        _ => Ok(None),
    }
}

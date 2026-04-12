use std::fmt;

use crate::error::{Error, Result};

/// A validated SRA accession identifier.
///
/// Accessions follow the pattern: prefix (2-3 letters) + 6-9 digits.
/// Prefixes: SRR/ERR/DRR (runs), SRS/ERS/DRS (samples), SRP/ERP/DRP (projects).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Accession {
    pub prefix: AccessionPrefix,
    pub number: String,
    raw: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AccessionPrefix {
    /// NCBI Sequence Read Archive
    Srr,
    /// EBI European Nucleotide Archive
    Err,
    /// DDBJ
    Drr,
}

impl fmt::Display for Accession {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.raw)
    }
}

impl fmt::Display for AccessionPrefix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Srr => write!(f, "SRR"),
            Self::Err => write!(f, "ERR"),
            Self::Drr => write!(f, "DRR"),
        }
    }
}

/// Parse and validate an SRA run accession string.
pub fn parse(input: &str) -> Result<Accession> {
    let input = input.trim();
    if input.len() < 9 {
        return Err(Error::InvalidAccession(format!(
            "'{input}' is too short for an SRA run accession"
        )));
    }

    let prefix_str = &input[..3];
    let prefix = match prefix_str.to_uppercase().as_str() {
        "SRR" => AccessionPrefix::Srr,
        "ERR" => AccessionPrefix::Err,
        "DRR" => AccessionPrefix::Drr,
        _ => {
            return Err(Error::InvalidAccession(format!(
                "'{input}' has unrecognized prefix '{prefix_str}' (expected SRR/ERR/DRR)"
            )));
        }
    };

    let number = &input[3..];
    if number.len() < 6 || number.len() > 9 {
        return Err(Error::InvalidAccession(format!(
            "'{input}' has {}-digit number (expected 6-9)",
            number.len()
        )));
    }
    if !number.chars().all(|c| c.is_ascii_digit()) {
        return Err(Error::InvalidAccession(format!(
            "'{input}' contains non-digit characters after prefix"
        )));
    }

    Ok(Accession {
        prefix,
        number: number.to_string(),
        raw: format!("{prefix}{number}"),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_accessions() {
        let acc = parse("SRR000001").unwrap();
        assert_eq!(acc.prefix, AccessionPrefix::Srr);
        assert_eq!(acc.number, "000001");
        assert_eq!(acc.to_string(), "SRR000001");

        let acc = parse("ERR1234567").unwrap();
        assert_eq!(acc.prefix, AccessionPrefix::Err);

        let acc = parse("DRR123456789").unwrap();
        assert_eq!(acc.prefix, AccessionPrefix::Drr);
        assert_eq!(acc.number, "123456789");
    }

    #[test]
    fn case_insensitive() {
        let acc = parse("srr000001").unwrap();
        assert_eq!(acc.prefix, AccessionPrefix::Srr);
    }

    #[test]
    fn trims_whitespace() {
        let acc = parse("  SRR000001  ").unwrap();
        assert_eq!(acc.to_string(), "SRR000001");
    }

    #[test]
    fn rejects_bad_prefix() {
        assert!(parse("XRR000001").is_err());
    }

    #[test]
    fn rejects_short_number() {
        assert!(parse("SRR12345").is_err());
    }

    #[test]
    fn rejects_non_digits() {
        assert!(parse("SRR00000a").is_err());
    }
}

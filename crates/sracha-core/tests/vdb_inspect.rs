//! Integration tests for the generic VDB inspector against a local fixture.
//!
//! Marked `#[ignore]` because the fixture lives outside the crate (in
//! `validation/`) and isn't packaged with `cargo publish`. Run with:
//!
//! ```bash
//! cargo test -p sracha-core -- --ignored vdb_inspect
//! ```

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sracha_core::vdb::inspect::{self, VdbKind};
use sracha_core::vdb::kar::KarArchive;
use sracha_core::vdb::metadata;

fn fixture(name: &str) -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation")
        .join(name);
    p.exists().then_some(p)
}

#[test]
#[ignore]
fn srr2584863_info_smoke() {
    let path = fixture("SRR2584863.sra").expect("validation/SRR2584863.sra");
    let f = File::open(&path).unwrap();
    let mut kar = KarArchive::open(BufReader::new(f)).unwrap();

    let kind = inspect::detect_kind(&kar).unwrap();
    assert!(matches!(kind, VdbKind::Database | VdbKind::Table));

    let cols = inspect::list_columns(&kar, None).unwrap();
    assert!(cols.iter().any(|c| c == "READ"), "READ column missing");

    let report = inspect::gather_info(&mut kar, &path).unwrap();
    assert!(report.primary_row_count().unwrap() > 0);

    if let Some(nodes) = inspect::read_table_metadata(&mut kar, None) {
        let _ = metadata::schema_text(&nodes);
        let _ = metadata::load_timestamp(&nodes);
    }
}

#[test]
#[ignore]
fn srr2584863_id_range_matches_row_count() {
    let path = fixture("SRR2584863.sra").expect("validation/SRR2584863.sra");
    let f = File::open(&path).unwrap();
    let mut kar = KarArchive::open(BufReader::new(f)).unwrap();
    let (first, count) = inspect::id_range(&mut kar, &path, None, None).unwrap();
    assert!(count > 0);
    assert_eq!(first, 1, "row IDs are 1-indexed");
}

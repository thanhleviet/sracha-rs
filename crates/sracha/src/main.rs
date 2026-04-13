mod cli;
mod style;

use anyhow::Result;
use clap::Parser;
use tracing_subscriber::EnvFilter;

use cli::{Cli, Command};
use sracha_core::sdl::ResolvedAccession;

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();

    // Set up tracing
    let filter = match (cli.quiet, cli.verbose) {
        (true, _) => "error",
        (_, 0) => "info",
        (_, 1) => "debug",
        (_, 2) => "trace",
        _ => "trace",
    };
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(filter)),
        )
        .with_writer(std::io::stderr)
        .init();

    match cli.command {
        Command::Fetch(args) => {
            tracing::info!("fetching {} accession(s)", args.accessions.len());
            for acc in &args.accessions {
                let acc = sracha_core::accession::parse(acc)?;
                tracing::info!("resolved accession: {acc}");
            }
            anyhow::bail!("fetch not yet implemented")
        }
        Command::Fastq(args) => {
            tracing::info!("converting {} input(s) to FASTQ", args.inputs.len());

            let split_mode = match args.split {
                cli::SplitMode::Split3 => sracha_core::fastq::SplitMode::Split3,
                cli::SplitMode::SplitFiles => sracha_core::fastq::SplitMode::SplitFiles,
                cli::SplitMode::SplitSpot => sracha_core::fastq::SplitMode::SplitSpot,
                cli::SplitMode::Interleaved => sracha_core::fastq::SplitMode::Interleaved,
            };

            for input in &args.inputs {
                let sra_path = std::path::Path::new(input);
                if !sra_path.exists() {
                    eprintln!("{} file not found: {}", style::error_label("error:"), style::path(input));
                    continue;
                }

                let pipeline_config = sracha_core::pipeline::PipelineConfig {
                    output_dir: args.output_dir.clone(),
                    split_mode,
                    gzip: !args.no_gzip,
                    gzip_level: args.gzip_level,
                    threads: args.threads,
                    connections: 1,
                    skip_technical: !args.include_technical,
                    min_read_len: args.min_read_len,
                    force: args.force,
                    progress: args.progress,
                };

                let stats = sracha_core::pipeline::run_fastq(sra_path, None, &pipeline_config)?;

                eprintln!(
                    "{}: {} spots, {} reads written",
                    style::header(&stats.accession),
                    style::count(stats.spots_read),
                    style::count(stats.reads_written),
                );
                for path in &stats.output_files {
                    eprintln!("  -> {}", style::path(path.display()));
                }
            }

            Ok(())
        }
        Command::Get(args) => {
            tracing::info!("get {} accession(s) -> FASTQ", args.accessions.len());

            let split_mode = match args.split {
                cli::SplitMode::Split3 => sracha_core::fastq::SplitMode::Split3,
                cli::SplitMode::SplitFiles => sracha_core::fastq::SplitMode::SplitFiles,
                cli::SplitMode::SplitSpot => sracha_core::fastq::SplitMode::SplitSpot,
                cli::SplitMode::Interleaved => sracha_core::fastq::SplitMode::Interleaved,
            };

            let sdl_client = sracha_core::sdl::SdlClient::new();

            for acc_str in &args.accessions {
                let acc = sracha_core::accession::parse(acc_str)?;
                let resolved = sdl_client.resolve_one(&acc.to_string()).await?;

                let pipeline_config = sracha_core::pipeline::PipelineConfig {
                    output_dir: args.output_dir.clone(),
                    split_mode,
                    gzip: !args.no_gzip,
                    gzip_level: args.gzip_level,
                    threads: args.threads,
                    connections: args.connections,
                    skip_technical: !args.include_technical,
                    min_read_len: args.min_read_len,
                    force: args.force,
                    progress: args.progress,
                };

                let stats = sracha_core::pipeline::run_get(&resolved, &pipeline_config).await?;

                let pct = if stats.total_sra_size > 0 {
                    stats.bytes_downloaded as f64 / stats.total_sra_size as f64 * 100.0
                } else {
                    100.0
                };
                eprintln!(
                    "{}: {} spots, {} reads written, {} downloaded ({} of {})",
                    style::header(&stats.accession),
                    style::count(stats.spots_read),
                    style::count(stats.reads_written),
                    style::value(format_size(stats.bytes_downloaded)),
                    style::percentage(format!("{pct:.1}%")),
                    style::value(format_size(stats.total_sra_size)),
                );
                for path in &stats.output_files {
                    eprintln!("  -> {}", style::path(path.display()));
                }
            }

            Ok(())
        }
        Command::Info(args) => {
            let client = sracha_core::sdl::SdlClient::new();
            for (i, acc_str) in args.accessions.iter().enumerate() {
                if i > 0 {
                    println!();
                }
                let acc = match sracha_core::accession::parse(acc_str) {
                    Ok(a) => a,
                    Err(e) => {
                        eprintln!("{} {acc_str}: {e}", style::error_label("error:"));
                        continue;
                    }
                };
                match client.resolve_one(&acc.to_string()).await {
                    Ok(resolved) => print_resolved(&resolved),
                    Err(e) => eprintln!("{} {acc_str}: {e}", style::error_label("error:")),
                }
            }
            Ok(())
        }
    }
}

fn print_resolved(resolved: &ResolvedAccession) {
    let f = &resolved.sra_file;

    println!("{}", style::header(&resolved.accession));
    println!("  {}    {}", style::label("Size:"), style::value(format_size(f.size)));
    println!(
        "  {}     {}",
        style::label("MD5:"),
        style::value(f.md5.as_deref().unwrap_or("not provided"))
    );
    println!(
        "  {}    {}",
        style::label("Lite:"),
        style::value(if f.is_lite { "yes" } else { "no" })
    );
    println!("  {} {}", style::label("Mirrors:"), style::count(f.mirrors.len()));
    for m in &f.mirrors {
        println!("    [{}] {}", style::value(&m.service), style::path(&m.url));
    }

    if let Some(ref vdb) = resolved.vdbcache_file {
        println!(
            "  {} yes ({}, {} mirrors)",
            style::label("VDBcache:"),
            style::value(format_size(vdb.size)),
            style::count(vdb.mirrors.len())
        );
    }
}

fn format_size(bytes: u64) -> String {
    const KIB: u64 = 1024;
    const MIB: u64 = 1024 * KIB;
    const GIB: u64 = 1024 * MIB;
    const TIB: u64 = 1024 * GIB;

    if bytes >= TIB {
        format!("{:.2} TiB", bytes as f64 / TIB as f64)
    } else if bytes >= GIB {
        format!("{:.2} GiB", bytes as f64 / GIB as f64)
    } else if bytes >= MIB {
        format!("{:.2} MiB", bytes as f64 / MIB as f64)
    } else if bytes >= KIB {
        format!("{:.2} KiB", bytes as f64 / KIB as f64)
    } else {
        format!("{bytes} B")
    }
}

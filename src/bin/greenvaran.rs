use clap::{Parser, Subcommand};
use greenvaran::cli::{LOGO, VERSION};
use greenvaran::querytab::{run_querytab, QuerytabArgs};
use greenvaran::smallvars::{run_smallvars, SmallvarsArgs};
use greenvaran::sv::{run_sv, SvArgs};

#[derive(Debug, Parser)]
#[command(name = "greenvaran", disable_version_flag = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Annotate small variants VCF
    Smallvars(SmallvarsArgs),
    /// Annotate structural variants VCF
    Sv(SvArgs),
    /// Generate a table for greendb_query from VCF
    Querytab(QuerytabArgs),
    /// Print version
    Version,
}

fn print_dispatcher_help(unknown: Option<&str>) {
    eprintln!("\nCommands: ");
    println!("  {:<13}:   {}", "smallvars", "Annotate small variants VCF");
    println!("  {:<13}:   {}", "sv", "Annotate structural variants VCF");
    println!(
        "  {:<13}:   {}",
        "querytab", "Generate a table for greendb_query from VCF"
    );
    println!("  {:<13}:   {}", "version", "Print version");
    if let Some(cmd) = unknown {
        println!("unknown program '{cmd}'");
    }
}

fn main() {
    let raw_args: Vec<String> = std::env::args().collect();
    let known = ["smallvars", "sv", "querytab", "version"];

    if raw_args.len() < 2 {
        print_dispatcher_help(None);
        return;
    }

    let first = raw_args[1].as_str();
    if !known.contains(&first) {
        if first == "-h" || first == "-help" || first == "--help" {
            print_dispatcher_help(None);
        } else {
            print_dispatcher_help(Some(first));
        }
        return;
    }

    println!("{LOGO}");

    // Match Nim behavior: bare subcommand with no options shows subcommand --help
    if raw_args.len() == 2 && first != "version" {
        let mut fake = raw_args.clone();
        fake.push("--help".to_string());
        let _ = Cli::try_parse_from(fake).map_err(|e| e.print());
        return;
    }

    let cli = Cli::parse();
    let res = match cli.command {
        Commands::Smallvars(args) => run_smallvars(args),
        Commands::Sv(args) => run_sv(args),
        Commands::Querytab(args) => run_querytab(args),
        Commands::Version => {
            println!("GREEN-VARAN version {VERSION}");
            Ok(())
        }
    };

    if let Err(err) = res {
        eprintln!("FATAL: {err:#}");
        std::process::exit(1);
    }
}

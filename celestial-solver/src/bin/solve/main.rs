mod cli;
mod debug_images;
mod report;
mod run;
mod timing;

use clap::Parser;
use std::process::ExitCode;

fn main() -> ExitCode {
    let cli = cli::Cli::parse();
    init_logging(cli.verbose);
    match run::run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        // The same text and exit status as returning the error from `main`.
        Err(e) => {
            eprintln!("Error: {e:?}");
            ExitCode::FAILURE
        }
    }
}

fn init_logging(verbose: bool) {
    let level = if verbose {
        log::LevelFilter::Debug
    } else {
        log::LevelFilter::Warn
    };
    env_logger::Builder::new()
        .filter_level(level)
        .parse_default_env()
        .format_target(false)
        .format_timestamp(None)
        .init();
}

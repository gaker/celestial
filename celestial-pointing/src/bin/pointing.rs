use celestial_pointing::commands::{self, CommandOutput};
use celestial_pointing::session::Session;
use rustyline::completion::{Completer, Pair};
use rustyline::error::ReadlineError;
use rustyline::highlight::Highlighter;
use rustyline::hint::Hinter;
use rustyline::validate::Validator;
use rustyline::{Helper, Editor};
use std::fs;
use std::path::{Path, PathBuf};

fn history_path() -> PathBuf {
    dirs::home_dir()
        .unwrap_or_default()
        .join(".celestial_pointing_history")
}

struct PointingHelper {
    commands: Vec<String>,
    terms: Vec<String>,
}

impl PointingHelper {
    fn new() -> Self {
        Self {
            commands: [
                "APPLY", "CORRECT", "INDAT", "INMOD", "OUTMOD", "USE",
                "LOSE", "FIT", "CLIST", "SLIST", "SHOW", "RESET", "MASK",
                "UNMASK", "MVET", "OUTL", "FIX", "UNFIX", "PARALLEL",
                "CHAIN", "ADJUST", "FAUTO", "OPTIMAL", "LST", "PREDICT",
                "GSCAT", "GDIST", "GMAP", "GHA", "GDEC", "GHYST",
                "HELP", "QUIT",
            ].iter().map(|s| s.to_string()).collect(),
            terms: [
                "IH", "ID", "CH", "NP", "MA", "ME", "TF", "TX", "DAF", "FO",
                "HCES", "HCEC", "DCES", "DCEC",
                "IA", "IE", "CA", "NPAE", "AN", "AW",
            ].iter().map(|s| s.to_string()).collect(),
        }
    }
}

fn split_path_prefix(partial: &str) -> (&Path, &str) {
    if partial.is_empty() {
        return (Path::new("."), "");
    }
    let path = Path::new(partial);
    if partial.ends_with(std::path::is_separator) {
        return (path, "");
    }
    match (path.parent(), path.file_name()) {
        (Some(p), Some(f)) => {
            let dir = if p.as_os_str().is_empty() { Path::new(".") } else { p };
            (dir, f.to_str().unwrap_or(""))
        }
        _ => (Path::new("."), partial),
    }
}

fn complete_path(partial: &str) -> Vec<Pair> {
    let (dir, prefix) = split_path_prefix(partial);
    let entries = match fs::read_dir(dir) {
        Ok(entries) => entries,
        Err(_) => return vec![],
    };
    let base = if partial.is_empty() || !partial.contains(std::path::is_separator) {
        String::new()
    } else if partial.ends_with(std::path::is_separator) {
        partial.to_string()
    } else {
        let i = partial.rfind(std::path::is_separator).map_or(0, |i| i + 1);
        partial[..i].to_string()
    };
    entries
        .filter_map(|e| e.ok())
        .filter_map(|e| build_path_pair(&e, prefix, &base))
        .collect()
}

fn build_path_pair(entry: &fs::DirEntry, prefix: &str, base: &str) -> Option<Pair> {
    let name = entry.file_name().into_string().ok()?;
    if !name.starts_with(prefix) {
        return None;
    }
    let suffix = if entry.path().is_dir() { std::path::MAIN_SEPARATOR_STR } else { "" };
    Some(Pair {
        display: format!("{}{}", name, suffix),
        replacement: format!("{}{}{}", base, name, suffix),
    })
}

impl Completer for PointingHelper {
    type Candidate = Pair;

    fn complete(
        &self,
        line: &str,
        pos: usize,
        _ctx: &rustyline::Context<'_>,
    ) -> rustyline::Result<(usize, Vec<Pair>)> {
        let up_to = &line[..pos];
        let words: Vec<&str> = up_to.split_whitespace().collect();

        if words.is_empty() || (words.len() == 1 && !up_to.ends_with(' ')) {
            let prefix = words.first().map_or("", |s| *s).to_uppercase();
            let start = up_to.rfind(char::is_whitespace).map_or(0, |i| i + 1);
            let matches: Vec<Pair> = self.commands.iter()
                .filter(|c| c.starts_with(&prefix))
                .map(|c| Pair { display: c.clone(), replacement: c.clone() })
                .collect();
            Ok((start, matches))
        } else {
            let cmd = words[0].to_uppercase();
            if matches!(cmd.as_str(), "INDAT" | "INMOD" | "OUTMOD" | "GSCAT" | "GDIST" | "GMAP" | "GHA" | "GDEC" | "GHYST") {
                let partial = if up_to.ends_with(' ') { "" } else { words.last().copied().unwrap_or("") };
                let start = up_to.rfind(char::is_whitespace).map_or(0, |i| i + 1);
                let matches = complete_path(partial);
                Ok((start, matches))
            } else if matches!(cmd.as_str(), "USE" | "LOSE" | "FIX" | "UNFIX" | "PARALLEL" | "CHAIN") {
                let prefix = words.last().map_or("", |s| *s).to_uppercase();
                let start = up_to.rfind(char::is_whitespace).map_or(0, |i| i + 1);
                let matches: Vec<Pair> = self.terms.iter()
                    .filter(|t| t.starts_with(&prefix))
                    .map(|t| Pair { display: t.clone(), replacement: t.clone() })
                    .collect();
                Ok((start, matches))
            } else {
                Ok((pos, vec![]))
            }
        }
    }
}

impl Hinter for PointingHelper {
    type Hint = String;
}
impl Highlighter for PointingHelper {}
impl Validator for PointingHelper {}
impl Helper for PointingHelper {}

fn main() {
    println!("celestial-pointing v{}", env!("CARGO_PKG_VERSION"));
    println!("Type HELP for commands, Ctrl-D to exit\n");

    let helper = PointingHelper::new();
    let mut rl = match Editor::with_config(
        rustyline::Config::builder()
            .auto_add_history(true)
            .build(),
    ) {
        Ok(rl) => rl,
        Err(e) => {
            eprintln!("Failed to initialize editor: {}", e);
            return;
        }
    };
    rl.set_helper(Some(helper));

    let history = history_path();
    let _ = rl.load_history(&history);

    let mut session = Session::new();

    loop {
        let prompt = ">> ".to_string();

        match rl.readline(&prompt) {
            Ok(line) => {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                if line.eq_ignore_ascii_case("QUIT") {
                    println!("Bye!");
                    break;
                }
                match commands::dispatch(&mut session, line) {
                    Ok(output) => print_output(output),
                    Err(e) => eprintln!("Error: {}", e),
                }
            }
            Err(ReadlineError::Interrupted) => {
                println!("^C");
                continue;
            }
            Err(ReadlineError::Eof) => {
                println!("Bye!");
                break;
            }
            Err(e) => {
                eprintln!("Error: {}", e);
                break;
            }
        }
    }

    let _ = rl.save_history(&history);
}

fn print_output(output: CommandOutput) {
    match output {
        CommandOutput::Text(s) => println!("{}", s),
        CommandOutput::Table { headers, rows } => print_table(&headers, &rows),
        CommandOutput::FitDisplay(fit) => print_fit(&fit),
        CommandOutput::None => {}
    }
}

fn print_fit(fit: &commands::FitDisplay) {
    println!("\n       coeff          value      sigma\n");
    for (i, name) in fit.term_names.iter().enumerate() {
        println!(
            "{:3}  {:>6}    {:>12.2}  {:>9.3}",
            i + 1, name, fit.coefficients[i], fit.sigma[i],
        );
    }
    println!("\nSky RMS = {:.2}\"\n", fit.sky_rms);
}

fn print_table(headers: &[String], rows: &[Vec<String>]) {
    let widths: Vec<usize> = (0..headers.len())
        .map(|i| {
            let hw = headers[i].len();
            let rw = rows.iter()
                .map(|r| r.get(i).map_or(0, |s| s.len()))
                .max()
                .unwrap_or(0);
            hw.max(rw)
        })
        .collect();

    for (i, h) in headers.iter().enumerate() {
        print!("{:>width$}  ", h, width = widths[i]);
    }
    println!();

    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            print!("{:>width$}  ", cell, width = widths[i]);
        }
        println!();
    }
}

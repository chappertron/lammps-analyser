/// TODO:
/// - [x] Scan all identifiers in one pass?
/// - [ ] Add incrental parsing
/// - [x] Make the output look a bit prettier
/// - [ ] Add file snippets to show where errors are
/// - [ ] Add a test suite
/// - [ ] Create an issue abstraction to handle all errors and warnings
/// - [ ] Extract core behaviour into a library
use anyhow::Result;
use ariadne::{Label, Report, ReportKind, Source};
use clap::Parser as ClapParser;
use lammps_analyser::{
    check_styles::check_styles, diagnostic_report::ReportSimple, error_finder::ErrorFinder,
    identifinder::IdentiFinder, lammps_errors::LammpsError,
};
use owo_colors::OwoColorize;
use std::{
    fs::File,
    io::{BufReader, Read},
};
use tree_sitter::Parser;

#[derive(Debug, clap::Parser)]

struct Cli {
    source: String,
    /// Output the parsed tree to a dot file
    #[clap(long)]
    output_tree: bool,
    /// Prints more complicated reports that highlight error locations
    #[clap(long)]
    output_reports: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let file = File::open(&cli.source)?;

    let mut file = BufReader::new(file);

    let mut source_code = String::new();

    file.read_to_string(&mut source_code)?;

    let source_bytes = source_code.as_bytes();

    let mut issues: Vec<LammpsError> = Vec::new();

    let mut parser = Parser::new();

    parser
        .set_language(tree_sitter_lammps::language())
        .expect("Could not load language");

    let tree = parser.parse(source_bytes, None).unwrap();

    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        tree.print_dot_graph(&dot_file);
    }

    let mut identifinder = IdentiFinder::new(&tree, source_bytes)?;

    let undefined_fixes = match identifinder.check_symbols() {
        Ok(()) => vec![],
        Err(v) => v,
    };

    let mut error_finder = ErrorFinder::new()?;
    _ = error_finder.find_syntax_errors(&tree, source_bytes)?;
    error_finder.find_missing_nodes(&tree)?;
    let syntax_errors = error_finder.syntax_errors();

    for err in syntax_errors {
        // ruff_test.py:3:1: F821 Undefined name `hello`
        println!("{}:{}", cli.source.bold(), err.make_simple_report());
    }

    // this doesn't work. Need to compare the names!

    for err in &undefined_fixes {
        // ruff_test.py:3:1: F821 Undefined name `hello`
        // println!("{}",std::str::from_utf8(source_code[ident.start_byte..ident.end_byte])?.underline());
        if cli.output_reports {
            Report::build(ReportKind::Error, &cli.source, err.ident.start_byte)
                .with_label(
                    Label::new((&cli.source, err.ident.start_byte..err.ident.end_byte))
                        .with_message(format!("{}", err)),
                )
                .finish()
                .print((
                    &cli.source,
                    Source::from(std::str::from_utf8(source_bytes)?),
                ))?;
        }
        println!("{}:{}", cli.source.bold(), err.make_simple_report());
    }

    let invalid_styles = check_styles(&tree, source_bytes)?;

    for err in &invalid_styles {
        println!("{}:{}", cli.source.bold(), err.make_simple_report());
    }
    // TODO Check if any warnings or errors are found!!!

    issues.extend(syntax_errors.iter().map(|x| x.clone().into()));
    issues.extend(undefined_fixes.iter().map(|x| x.clone().into()));
    issues.extend(invalid_styles.iter().map(|x| x.clone().into()));

    if !issues.is_empty() {
        let n_errors = issues.len();
        println!(
            "{}: {} error{} found 😞",
            cli.source.bold(),
            n_errors.bright_red(),
            if n_errors != 1 { "s" } else { "" },
        );
        Ok(())
    } else {
        println!("All Good 😊");
        Ok(())
    }
}

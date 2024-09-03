/// TODO:
/// - [ ] Add incrental parsing
/// - [ ] Add file snippets to show where errors are.
/// - [ ] Add a test suite
/// - [ ] Create an issue abstraction to handle all errors and warnings
/// - [ ] Extract core behaviour into a library
use anyhow::{Context, Result};

use clap::Parser as ClapParser;
use lammps_analyser::{diagnostic_report::ReportSimple, input_script};
use owo_colors::OwoColorize;
use std::fs::File;
use tree_sitter::Parser;

#[derive(Debug, clap::Parser)]
struct Cli {
    /// LAMMPS input script to check
    source: String,
    /// Output the parsed tree to a dot file.
    /// This can be used to visualise the tree with e.g. graphviz
    #[clap(long)]
    output_tree: bool,
    /// Prints more complicated reports that highlight error locations
    #[clap(long)]
    output_reports: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let source_code = std::fs::read_to_string(&cli.source)?;
    let mut parser = Parser::new();

    parser
        .set_language(tree_sitter_lammps::language())
        .context("Could not load tree-sitter language")?;

    let state = input_script::InputScript::new(&source_code, &mut parser)?;

    // Output a syntax tree for debugging.
    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        state.tree.print_dot_graph(&dot_file);
    }

    // Errors
    if !state.ast_errors.is_empty() {
        for error in &state.ast_errors {
            println!("{}:{}", cli.source.bold(), error.make_simple_report());
        }
    }

    for issue in &state.issues {
        println!("{}", issue.make_simple_report());
    }

    for diagnostic in &state.diagnostics {
        println!("{}", diagnostic.make_simple_report());
    }
    if !state.issues.is_empty() && !state.diagnostics.is_empty() {
        // TODO:   Count warnings separately!!!
        let n_errors = state.issues.len() + state.diagnostics.len();
        println!(
            "{}: {} issues{} found 😞",
            cli.source.bold(),
            n_errors.bright_red(),
            if n_errors == 1 { "" } else { "s" },
        );
    } else {
        println!("All Good 😊");
    }
    Ok(())
}

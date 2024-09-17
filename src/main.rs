/// TODO:
/// - [ ] Add incrental parsing
/// - [ ] Add file snippets to show where errors are.
/// - [ ] Add a test suite
/// - [ ] Create an issue abstraction to handle all errors and warnings
/// - [ ] Extract core behaviour into a library
use anyhow::{Context, Result};

use clap::Parser as ClapParser;
use lammps_analyser::{diagnostic_report::FileNameReport, input_script};
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

    let state = input_script::InputScript::new(&source_code)?;

    // Output a syntax tree for debugging.
    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        state.tree.print_dot_graph(&dot_file);
    }

    for diagnostic in &state.diagnostics {
        println!("{}", diagnostic.make_file_name_report(&cli.source));
    }
    if !state.diagnostics.is_empty() {
        // TODO:   Count warnings separately!!!
        let n_errors = state.diagnostics.len();
        println!(
            "{}: {} issue{} found 😞",
            cli.source.bold(),
            n_errors.bright_red(),
            if n_errors == 1 { "" } else { "s" },
        );

        std::process::exit(72);
    } else {
        println!("All Good 😊");
        std::process::exit(0);
    }
}

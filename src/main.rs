/// TODO:
/// - [x] Scan all identifiers in one pass?
/// - [ ] Add incrental parsing
/// - [x] Make the output look a bit prettier
/// - [ ] Add file snippets to show where errors are
/// - [ ] Add a test suite
/// - [ ] Create an issue abstraction to handle all errors and warnings
/// - [ ] Extract core behaviour into a library
use anyhow::Result;

use clap::Parser as ClapParser;
use lammps_analyser::{
    ast::{ts_to_ast, CommandType, NamedCommand},
    check_commands,
    check_styles::check_styles,
    diagnostic_report::ReportSimple,
    error_finder::ErrorFinder,
    identifinder::{unused_variables, IdentiFinder},
    issues::Issue,
    lammps_errors::{LammpsError, Warnings},
};
use owo_colors::OwoColorize;
use std::{
    fs::File,
    io::{BufReader, Read},
};
use tree_sitter::Parser;

#[derive(Debug, clap::Parser)]

struct Cli {
    /// LAMMPS input script to check
    source: String,
    /// Output the parsed tree to a dot file
    /// This can be used to visualise the tree with e.g. graphviz
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

    let mut issues: Vec<Issue> = Vec::new();

    let mut parser = Parser::new();

    parser
        .set_language(tree_sitter_lammps::language())
        .expect("Could not load language");

    let tree = parser.parse(source_bytes, None).unwrap();

    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        tree.print_dot_graph(&dot_file);
    }

    let ast = ts_to_ast(&tree, source_bytes);

    // Somewhat gracefully exit
    if let Err(e) = &ast {
        println!("{}:{}", cli.source.bold(), e.make_simple_report());
    }
    let ast = ast?;

    // Parsing fixes

    let parsed_fixes = ast
        .commands
        .iter()
        .filter_map(|command| {
            if let CommandType::NamedCommand(NamedCommand::Fix(fix)) = &command.command_type {
                Some(check_commands::check_fix::check_fix(fix))
            } else {
                None
            }
        })
        .filter_map(|x| x.err())
        .collect::<Vec<_>>();

    let identifinder = IdentiFinder::new(&tree, source_bytes)?;

    let undefined_fixes = match identifinder.check_symbols() {
        Ok(()) => vec![],
        Err(v) => v,
    };

    let mut error_finder = ErrorFinder::new()?;
    _ = error_finder.find_syntax_errors(&tree, source_bytes)?;
    error_finder.find_missing_nodes(&tree)?;
    let syntax_errors = error_finder.syntax_errors();

    // TODO: Check if any warnings or errors are found!!!

    let invalid_styles = check_styles(&tree, source_bytes)?;
    issues.extend(
        syntax_errors
            .iter()
            .map(|x| LammpsError::from(x.clone()).into()),
    );
    issues.extend(
        undefined_fixes
            .into_iter()
            .map(|x| LammpsError::from(x.clone()).into()),
    );
    issues.extend(
        invalid_styles
            .into_iter()
            .map(|x| LammpsError::from(x.clone()).into()),
    );
    issues.extend(
        unused_variables(identifinder.symbols())
            .into_iter()
            .map(|x| Warnings::from(x).into()),
    );

    issues.extend(
        parsed_fixes
            .into_iter()
            .map(|x| LammpsError::from(x).into()),
    );

    for issue in &issues {
        println!("{}", issue.make_simple_report());
    }
    if !issues.is_empty() {
        // TODO:  Don't count warnings as errors!!!
        let n_errors = issues.len();
        println!(
            "{}: {} error{} found 😞",
            cli.source.bold(),
            n_errors.bright_red(),
            if n_errors == 1 { "" } else { "s" },
        );
    } else {
        println!("All Good 😊");
    }
    Ok(())
}

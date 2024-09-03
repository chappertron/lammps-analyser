/// TODO:
/// - [ ] Add incrental parsing
/// - [ ] Add file snippets to show where errors are.
/// - [ ] Add a test suite
/// - [ ] Create an issue abstraction to handle all errors and warnings
/// - [ ] Extract core behaviour into a library
use anyhow::{Context, Result};

use clap::Parser as ClapParser;
use lammps_analyser::issues::Issue as ScriptIssue;
use lammps_analyser::{
    ast::{from_node::FromNodeError, ts_to_ast, Ast, CommandType, NamedCommand, PartialAst},
    check_commands,
    check_styles::check_styles,
    diagnostic_report::ReportSimple,
    diagnostics::{Diagnostic, Issue},
    error_finder::ErrorFinder,
    identifinder::{unused_variables, IdentiFinder},
    lammps_errors::{LammpsError, Warnings},
    spanned_error::SpannedError,
};
use owo_colors::OwoColorize;
use std::fs::File;
use tree_sitter::{Parser, Tree};

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

#[derive(Debug)] // TODO: Allow for implementing clone. Can't yet because of Query in Identifinder.
pub struct InputScript<'src> {
    pub source_code: &'src str,
    pub diagnostics: Vec<Diagnostic>,
    pub issues: Vec<ScriptIssue>, // TODO: Use the new issue trait?
    pub tree: Tree,
    pub ast: Ast,
    pub ast_errors: Option<Vec<SpannedError<FromNodeError>>>,
    pub identifinder: IdentiFinder,
    pub error_finder: ErrorFinder,
}

impl<'src> InputScript<'src> {
    /// Monolithic method that reads the lammps source code.
    /// Parser is taken as input rather than stored because it does not implement debug.
    fn new(source_code: &'src str, parser: &mut Parser) -> Result<Self> {
        let tree = parser
            .parse(source_code, None)
            .context("Failed to load the TS grammar.")?;

        let mut issues: Vec<ScriptIssue> = Vec::new();
        let mut diagnostics: Vec<Diagnostic> = Vec::new();

        let ast = ts_to_ast(&tree, source_code);
        // Somewhat gracefully exit

        let (ast, ast_errors) = match ast {
            Ok(ast) => (ast, None),
            Err(PartialAst { ast, errors }) => (ast, Some(errors)),
        };

        // Checking fix arguments
        let fix_errors = ast
            .commands
            .iter()
            .filter_map(|command| {
                // TODO: Use a checkcommand function that checks all command types.
                if let CommandType::NamedCommand(NamedCommand::Fix(fix)) = &command.command_type {
                    Some(check_commands::fixes::check_fix(fix))
                } else if let CommandType::NamedCommand(NamedCommand::Compute(compute)) =
                    &command.command_type
                {
                    Some(check_commands::computes::check_compute(compute))
                } else {
                    None
                }
            })
            .filter_map(|x| x.err())
            .map(|issue| issue.diagnostic())
            .collect::<Vec<_>>();

        let identifinder = IdentiFinder::new(&tree, source_code)?;

        let undefined_fixes = match identifinder.check_symbols() {
            Ok(()) => vec![],
            Err(v) => v,
        };

        let mut error_finder = ErrorFinder::new()?;
        _ = error_finder.find_syntax_errors(&tree, source_code)?;
        error_finder.find_missing_nodes(&tree)?;
        let syntax_errors = error_finder.syntax_errors();

        let invalid_styles = check_styles(&tree, source_code)?;
        issues.extend(
            syntax_errors
                .iter()
                .map(|x| LammpsError::from(x.clone()).into()),
        );
        issues.extend(
            undefined_fixes
                .into_iter()
                .map(|x| LammpsError::from(x).into()),
        );
        issues.extend(
            invalid_styles
                .into_iter()
                .map(|x| LammpsError::from(x).into()),
        );
        issues.extend(
            unused_variables(identifinder.symbols())
                .into_iter()
                .map(|x| Warnings::from(x).into()),
        );

        diagnostics.extend(fix_errors);

        Ok(Self {
            source_code,
            tree,
            issues,
            ast,
            ast_errors,
            diagnostics,
            identifinder,
            error_finder,
        })
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let source_code = std::fs::read_to_string(&cli.source)?;
    let mut parser = Parser::new();

    parser
        .set_language(tree_sitter_lammps::language())
        .context("Could not load tree-sitter language")?;

    let state = InputScript::new(&source_code, &mut parser)?;

    // Output a syntax tree for debugging.
    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        state.tree.print_dot_graph(&dot_file);
    }

    // Errors
    if let Some(errors) = state.ast_errors {
        for error in errors {
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

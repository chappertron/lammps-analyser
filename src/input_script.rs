use std::fmt::Debug;

use crate::check_commands::{self};
use crate::error_finder::SyntaxError;
use crate::issues::Issue as ScriptIssue;
use crate::lammps_errors::Warnings;
use crate::utils;
use crate::{
    ast::from_node::FromNodeError, check_styles::check_styles, diagnostics::Issue,
    error_finder::ErrorFinder, spanned_error::SpannedError,
};

use crate::identifinder::unused_variables;

use crate::lammps_errors::LammpsError;

use crate::ast::{CommandType, PartialAst};

use crate::ast::ts_to_ast;

use anyhow::Result;

use crate::identifinder::IdentiFinder;

use crate::ast::Ast;

use crate::diagnostics::Diagnostic;
use anyhow::Context;
use tree_sitter::{Parser, Tree};

#[derive(Debug)] // TODO: Allow for implementing clone. Can't yet because of Query in Identifinder.
pub struct InputScript<'src> {
    // TODO: add a file name.
    pub source_code: &'src str,
    pub diagnostics: Vec<Diagnostic>,
    pub issues: Vec<ScriptIssue>, // TODO: Use the new issue trait?
    /// The tree_sitter concrete tree
    pub tree: Tree,
    /// The higher level abstract syntax tree
    pub ast: Ast,
    pub identifinder: IdentiFinder,
    pub error_finder: ErrorFinder,
    parser: LmpParser,
}

/// Implemented so `Debug` can be derived for `InputScript`
struct LmpParser(Parser);
impl Debug for LmpParser {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("LmpParser").field(&"Parser").finish()
    }
}

impl<'src> InputScript<'src> {
    /// Monolithic method that reads the lammps source code.
    /// Parser is taken as input rather than stored because it does not implement debug.
    pub fn new(source_code: &'src str) -> Result<Self> {
        let mut parser = utils::parsing::setup_parser();
        let tree = parser
            .parse(source_code, None)
            .context("Failed to load the TS grammar.")?;

        let mut issues: Vec<ScriptIssue> = Vec::new();
        let mut diagnostics: Vec<Diagnostic> = Vec::new();

        let ast = ts_to_ast(&tree, source_code);
        // Somewhat gracefully exit

        let (ast, ast_errors) = match ast {
            Ok(ast) => (ast, vec![]),
            Err(PartialAst { ast, errors }) => (ast, errors),
        };

        // Checking fix arguments
        let fix_errors = ast
            .commands
            .iter()
            .filter_map(|command| {
                // TODO: Use a checkcommand function that checks all command types.
                if let CommandType::Fix(fix) = &command.command_type {
                    Some(check_commands::fixes::check_fix(fix))
                } else if let CommandType::Compute(compute) = &command.command_type {
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

        diagnostics.extend(
            unused_variables(identifinder.symbols())
                .into_iter()
                .map(|x| Warnings::from(x).diagnostic()),
        );

        // TODO: convert to diagnostics
        diagnostics.extend(syntax_errors.iter().cloned().map(|x| x.diagnostic()));

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
            ast_errors
                .into_iter()
                .map(|x| LammpsError::from(SyntaxError::from(x)).into()),
        );

        diagnostics.extend(fix_errors);

        Ok(Self {
            parser: LmpParser(parser),
            source_code,
            tree,
            issues,
            ast,
            diagnostics,
            identifinder,
            error_finder,
        })
    }

    fn issues(&self) -> impl Iterator<Item = Diagnostic> {
        todo!();
        vec![].into_iter()
    }
}

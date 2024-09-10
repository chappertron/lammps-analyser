use std::fmt::Debug;

use crate::check_commands::{self};
use crate::error_finder::SyntaxError;
use crate::lammps_errors::Warnings;
use crate::lints::redefined_identifiers;
use crate::utils;
use crate::{check_styles::check_styles, diagnostics::Issue, error_finder::ErrorFinder};

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
    /// The tree_sitter concrete tree
    pub tree: Tree,
    /// The higher level abstract syntax tree
    pub ast: Ast,
    pub identifinder: IdentiFinder,
    pub error_finder: ErrorFinder,
    _parser: LmpParser,
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
                // TODO: Use a check command function that checks all command types.
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

        error_finder.find_syntax_errors(&tree, source_code)?;
        error_finder.find_missing_nodes(&tree)?;

        let syntax_errors = error_finder.syntax_errors();

        let invalid_styles = check_styles(&ast, &tree, source_code)?;

        // These first because they are likely least severe.
        diagnostics.extend(fix_errors);

        diagnostics.extend(ast.check_commands().map(|e| e.diagnostic()));

        diagnostics.extend(
            unused_variables(identifinder.symbols())
                .into_iter()
                .map(|x| Warnings::from(x).diagnostic()),
        );

        // TODO: return owned syntax errors to remove cloning
        diagnostics.extend(syntax_errors.iter().cloned().map(|x| x.diagnostic()));

        // TODO: convert to diagnostics
        diagnostics.extend(
            undefined_fixes
                .into_iter()
                .map(|x| LammpsError::from(x).diagnostic()),
        );

        diagnostics.extend(
            invalid_styles
                .into_iter()
                .map(|x| LammpsError::from(x).diagnostic()),
        );

        diagnostics.extend(
            ast_errors
                .into_iter()
                .map(|x| LammpsError::from(SyntaxError::from(x)).diagnostic()),
        );

        let mut script = Self {
            _parser: LmpParser(parser),
            source_code,
            tree,
            ast,
            diagnostics,
            identifinder,
            error_finder,
        };
        script.run_lints();
        Ok(script)
    }

    pub fn diagnostics(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics.iter()
    }

    pub fn run_lints(&mut self) {
        self.diagnostics.extend(
            redefined_identifiers(&self.ast, &self.identifinder.symbols())
                .map(|id| id.diagnostic()),
        )
    }
}

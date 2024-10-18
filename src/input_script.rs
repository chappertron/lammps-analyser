use std::fmt::Debug;

use crate::lints::redefined_identifiers;
use crate::utils;
use crate::{check_styles::check_styles, diagnostics::Issue, error_finder::ErrorFinder};

use crate::identifinder::unused_references;

use crate::ast::PartialAst;

use crate::ast::ts_to_ast;

use anyhow::Result;

use crate::identifinder::IdentiFinder;

use crate::ast::Ast;

use crate::diagnostics::Diagnostic;
use anyhow::Context;
use tree_sitter::{Parser, Tree};

#[derive(Debug)] // TODO: Allow for implementing clone. Can't yet because of QueryCursor in Identifinder.
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
    parser: LmpParser,
}

/// Implemented so `Debug` can be derived for `InputScript`
#[allow(dead_code)]
struct LmpParser(Parser);

impl Debug for LmpParser {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("LmpParser").field(&"Parser").finish()
    }
}

impl<'src> InputScript<'src> {
    /// Method that reads the LAMMPS source code and reports and stores errors.
    pub fn new(source_code: &'src str) -> Result<Self> {
        let mut input_script = Self::new_minimal(source_code)?;

        input_script.check();

        Ok(input_script)
    }

    /// Creates a new input script and simply parses. Does not run any checks.
    fn new_minimal(source_code: &'src str) -> Result<Self> {
        let mut parser = utils::parsing::setup_parser();

        let tree = parser
            .parse(source_code, None)
            .context("Failed to load the TS grammar.")?;

        let diagnostics: Vec<Diagnostic> = Vec::new();

        Ok(Self {
            source_code,
            diagnostics,
            tree,
            ast: Ast::default(),
            identifinder: IdentiFinder::new_no_parse(),
            error_finder: ErrorFinder::new(),
            parser: LmpParser(parser),
        })
    }

    pub fn reparse_and_check(&mut self, src: &str) {
        self.reparse(src).expect("failed to parse");
        self.check();
    }

    /// Re-obtain the tree-sitter tree
    fn reparse(&mut self, src: &str) -> Result<()> {
        let new_tree = self
            .parser
            .0
            .parse(src, Some(&self.tree))
            .context("failed to re-parse text")?;
        self.tree = new_tree;

        Ok(())
    }

    /// Reparse the AST and re-run the checks on the script
    fn check(&mut self) {
        self.diagnostics.clear();

        let (ast, ast_errors) = get_ast(&self.tree, self.source_code);

        self.ast = ast;
        self.diagnostics.extend(ast_errors);

        self.check_symbols();

        self.check_commands();

        self.check_styles();

        self.unused_references();

        self.find_syntax_errors();
        self.run_lints();
    }

    /// Iterator over all diagnostics that have been found in the file.
    pub fn diagnostics(&self) -> impl Iterator<Item = &Diagnostic> {
        self.diagnostics.iter()
    }

    fn run_lints(&mut self) {
        self.diagnostics.extend(
            redefined_identifiers(&self.ast, self.identifinder.symbols()).map(|id| id.diagnostic()),
        )
    }

    /// Helper method to extend the diagnostics
    fn extend(&mut self, iter: impl IntoIterator<Item = impl Issue>) {
        self.diagnostics
            .extend(iter.into_iter().map(|e| e.diagnostic()))
    }

    fn find_syntax_errors(&mut self) {
        let syntax_errors = self.error_finder.find_errors(&self.tree, self.source_code);
        self.diagnostics
            .extend(syntax_errors.iter().map(|x| x.diagnostic()));
    }

    fn check_symbols(&mut self) {
        match self.identifinder.find_symbols(&self.tree, self.source_code) {
            Ok(_) => (),
            Err(error) => {
                self.diagnostics.push(error.diagnostic());
                return;
            }
        }

        let undefined_symbols = match self.identifinder.check_symbols() {
            Ok(()) => vec![],
            Err(v) => v,
        };

        self.extend(undefined_symbols);
    }

    fn check_commands(&mut self) {
        // TODO: Can't use helper method due to lifetime issues
        self.diagnostics
            .extend(self.ast.check_commands().map(|e| e.diagnostic()));
    }

    fn unused_references(&mut self) {
        let unused = unused_references(self.identifinder.symbols());
        self.extend(unused);
    }

    fn check_styles(&mut self) {
        self.diagnostics.extend(
            check_styles(&self.ast, &self.tree, self.source_code)
                .into_iter()
                .map(|x| x.diagnostic()),
        );
    }
}

fn get_ast(tree: &Tree, source_code: &str) -> (Ast, impl Iterator<Item = Diagnostic>) {
    let ast = ts_to_ast(tree, source_code);

    let (ast, ast_errors) = match ast {
        Ok(ast) => (ast, vec![]),
        Err(PartialAst { ast, errors }) => (ast, errors),
    };
    (ast, ast_errors.into_iter().map(|x| x.diagnostic()))
}

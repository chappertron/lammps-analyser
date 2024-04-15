use crate::spans::Point;
use crate::{diagnostics::Issue, spans::Span};
use anyhow::{Context, Result};
use lsp_types::Diagnostic as LspDiagnostic;
use lsp_types::DiagnosticSeverity;
use owo_colors::OwoColorize;
use std::fmt::Debug;
use thiserror::Error;
use tree_sitter::{Query, QueryCursor, Tree, TreeCursor};

use crate::diagnostic_report::ReportSimple;

/// Finds syntax issues in the file.
pub struct ErrorFinder {
    pub query: Query,
    cursor: QueryCursor,
    pub syntax_errors: Vec<SyntaxError>,
}

impl ErrorFinder {
    pub fn new() -> Result<Self> {
        let query = Query::new(
            tree_sitter_lammps::language(),
            "
            (ERROR) @syntax_error
            ;(MISSING) @syntax_error
            ",
        )?;
        let cursor = QueryCursor::new();
        let syntax_errors = vec![];

        Ok(Self {
            query,
            cursor,
            syntax_errors,
        })
    }

    // The slice of found syntax errors.
    pub fn syntax_errors(&self) -> &[SyntaxError] {
        self.syntax_errors.as_ref()
    }

    /// Find any syntax errors in the provided tree-sitter syntax tree.
    pub fn find_syntax_errors(
        &mut self,
        tree: &Tree,
        source_code: impl AsRef<[u8]>,
    ) -> Result<&[SyntaxError]> {
        let source_code = source_code.as_ref();
        let matches = self
            .cursor
            .matches(&self.query, tree.root_node(), source_code);
        for mat in matches {
            let text = mat.captures[0]
                .node
                .utf8_text(source_code)
                .context("Failed to parse to valid utf-8")?;
            let start = mat.captures[0].node.start_position();
            let end = mat.captures[0].node.end_position();
            self.syntax_errors.push(SyntaxError::ParseError(ParseError {
                text: text.into(),
                start: start.into(),
                end: end.into(),
            }));
        }
        Ok(self.syntax_errors())
    }

    /// Tree-sitter can't currently query for missing nodes, so recursivley walking the tree instead
    /// Missing nodes are also not reported as errors, so this is needed.
    /// TODO: Walk back from missing nodes to work out the expected node?
    pub fn find_missing_nodes(&mut self, tree: &Tree) -> Result<&Vec<SyntaxError>> {
        fn recur_missing(cursor: &mut TreeCursor, missing_nodes: &mut Vec<Point>) {
            if cursor.node().child_count() == 0 {
                if cursor.node().is_missing() {
                    let node = cursor.node();
                    missing_nodes.push(node.start_position().into());
                }
            } else {
                // Go to the first child, then recur
                cursor.goto_first_child();
                recur_missing(cursor, missing_nodes);
                // Go to the  siblings
                while cursor.goto_next_sibling() {
                    recur_missing(cursor, missing_nodes);
                }
                cursor.goto_parent();
            }
        }

        let mut cursor = tree.root_node().walk();
        let mut missing_nodes = vec![];

        recur_missing(&mut cursor, &mut missing_nodes);
        self.syntax_errors.extend(
            missing_nodes.iter().map(|start_position| {
                SyntaxError::MissingToken(MissingToken::new(*start_position))
            }),
        );

        Ok(&self.syntax_errors)
    }
}

impl Debug for ErrorFinder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ErrorFinder")
            .field("syntax_errors", &self.syntax_errors)
            .finish()
    }
}

#[derive(Error, Debug, PartialEq, Eq, Clone, Hash)]
#[error("{0}")]
pub enum SyntaxError {
    /// A token was expected, but none was found.
    MissingToken(MissingToken),
    /// Some other parsing error.
    ParseError(ParseError),
}
impl ReportSimple for SyntaxError {
    fn make_simple_report(&self) -> String {
        match self {
            Self::ParseError(parse_error) => parse_error.make_simple_report(),
            Self::MissingToken(missing_token) => missing_token.make_simple_report(),
        }
    }
}

impl Issue for SyntaxError {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        match self {
            Self::MissingToken(err) => err.diagnostic(),
            Self::ParseError(err) => err.diagnostic(),
        }
    }
}

impl From<SyntaxError> for lsp_types::Diagnostic {
    fn from(value: SyntaxError) -> Self {
        match value {
            SyntaxError::ParseError(parse_error) => LspDiagnostic::new(
                lsp_types::Range {
                    start: parse_error.start.into_lsp_type(),
                    end: parse_error.end.into_lsp_type(),
                },
                Some(DiagnosticSeverity::ERROR),
                None,
                None,
                parse_error.to_string(),
                None,
                None,
            ),
            SyntaxError::MissingToken(missing_token) => LspDiagnostic::new(
                lsp_types::Range {
                    start: missing_token.start.into_lsp_type(),
                    end: missing_token.start.into_lsp_type(),
                },
                Some(DiagnosticSeverity::ERROR),
                None,
                None,
                missing_token.to_string(),
                None,
                None,
            ),
        }
    }
}

// Perhaps store message into this
#[derive(Debug, PartialEq, Eq, Clone, Hash, Error)]
#[error(" {} `{}`", "Invalid Syntax", text)]
pub struct ParseError {
    pub text: String,
    pub start: Point,
    pub end: Point,
}
impl ReportSimple for ParseError {
    fn make_simple_report(&self) -> String {
        format!(
            "{}:{}: {} `{}`",
            self.start.row + 1,
            self.start.column + 1,
            "Invalid Syntax:".bright_red(),
            self.text
        )
    }
}

impl Issue for ParseError {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        crate::diagnostics::Diagnostic {
            name: "Syntax error:".to_owned(),
            severity: crate::diagnostics::Severity::Warning,
            span: Span {
                start: self.start,
                end: self.end,
            },
            message: format!("Invalid Syntax `{}`", self.text),
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash, Error)]
#[error("Invalid Syntax: Missing Token")]
pub struct MissingToken {
    pub start: Point,
}

impl MissingToken {
    fn new(start: Point) -> MissingToken {
        MissingToken { start }
    }
}

impl Issue for MissingToken {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        crate::diagnostics::Diagnostic {
            name: "Missing Token".into(),
            severity: crate::diagnostics::Severity::Error,
            span: Span {
                start: self.start,
                end: self.start,
            },
            message: self.to_string(),
        }
    }
}

// TODO: remove this implementation. Only have it implemented on diagnostic????
impl ReportSimple for MissingToken {
    fn make_simple_report(&self) -> String {
        format!(
            "{}:{}: {} `{}`",
            self.start.row + 1,
            self.start.column + 1,
            "Invalid Syntax:".bright_red(),
            "Missing Token",
        )
    }
}

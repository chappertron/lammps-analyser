use anyhow::Result;
use lsp_types::{Diagnostic, DiagnosticSeverity, Position};
use owo_colors::{OwoColorize, Stream::Stdout};
use std::fmt::Debug;
use thiserror::Error;
use tree_sitter::{Point, Query, QueryCursor, Tree, TreeCursor};

pub struct ErrorFinder {
    pub query: Query,
    cursor: QueryCursor,
    syntax_errors: Vec<SyntaxError>,
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

    pub fn syntax_errors(&self) -> &[SyntaxError] {
        self.syntax_errors.as_ref()
    }

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
            let text = mat.captures[0].node.utf8_text(source_code)?;
            let start = mat.captures[0].node.start_position();
            let end = mat.captures[0].node.end_position();
            self.syntax_errors.push(SyntaxError::ParseError(ParseError {
                text: text.into(),
                start,
                end,
            }));
        }
        Ok(self.syntax_errors())
    }

    /// Tree-sitter can't currently query for missing nodes, so recursivley walking the tree instead
    /// Missing nodes are also not reported as errors, so this is needed?
    /// TODO Walk back from missing nodes to work out the proper node?
    pub fn find_missing_nodes(&mut self, tree: &Tree) -> Result<&Vec<SyntaxError>> {
        let mut cursor = tree.root_node().walk();
        let mut missing_nodes = vec![];
        fn recur_missing(cursor: &mut TreeCursor, missing_nodes: &mut Vec<Point>) {
            if cursor.node().child_count() == 0 {
                if cursor.node().is_missing() {
                    // println!("{} {}:{}","Missing Node:".red(),cursor.node().start_position().row+1,cursor.node().start_position().column+1);
                    let node = cursor.node();
                    missing_nodes.push(node.start_position());
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

        recur_missing(&mut cursor, &mut missing_nodes);
        self.syntax_errors
            .extend(missing_nodes.iter().map(|start_position| {
                SyntaxError::MissingToken(MissingToken {
                    start: *start_position,
                })
            }));

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
    MissingToken(MissingToken),
    ParseError(ParseError),
}
fn point_to_position(point: Point) -> Position {
    Position {
        line: point.row as u32,
        character: point.column as u32,
    }
}

impl From<SyntaxError> for Diagnostic {
    fn from(value: SyntaxError) -> Self {
        match value {
            SyntaxError::ParseError(parse_error) => Diagnostic::new(
                lsp_types::Range {
                    start: point_to_position(parse_error.start),
                    end: point_to_position(parse_error.end),
                },
                Some(DiagnosticSeverity::ERROR),
                None,
                None,
                parse_error.to_string(),
                None,
                None,
            ),
            SyntaxError::MissingToken(missing_token) => Diagnostic::new(
                lsp_types::Range {
                    start: point_to_position(missing_token.start),
                    end: point_to_position(missing_token.start),
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

// TODO: Work out how to get rid of the colour tags when not in a terminal
// Perhaps store message into this
#[derive(Debug, PartialEq, Eq, Clone, Hash, Error)]
#[error("{}:{}: {} `{}`",start.row+1,start.column+1,"Invalid Syntax".if_supports_color(Stdout,|text|text.bright_red()),text)]
pub struct ParseError {
    pub text: String,
    pub start: Point,
    pub end: Point,
}

#[derive(Debug, PartialEq, Eq, Clone, Hash, Error)]
#[error("{}:{}: {} {}",start.row+1,start.column+1,"Invalid Syntax".if_supports_color(Stdout,|text|text.bright_red()),"Missing Token")]
pub struct MissingToken {
    pub start: Point,
}

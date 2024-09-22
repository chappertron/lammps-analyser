use tree_sitter::{Node, Point, TreeCursor};

use crate::{diagnostics, spanned_error::SpannedError, utils::into_error::IntoError};
use thiserror::Error;

use super::expressions::ParseExprError;

/// Converts from a tree-sitter node to the type this trait is implemented for.
///
/// Requires providing the text source, which is needed for extracting names and identifiers.
pub trait FromNode: Sized {
    type Error;
    fn from_node(node: &Node, text: &str) -> Result<Self, Self::Error>;
}

// TODO: Just wrap expression parse error??
#[derive(Debug, Error, Clone, PartialEq, Eq)]
/// An error associated with converting a tree-sitter `Node`.
pub enum FromNodeError {
    #[error("Could not parse text as int {0}")]
    ParseIntError(std::num::ParseIntError),
    #[error("Could not parse text as float {0}")]
    ParseFloatError(std::num::ParseFloatError),
    #[error("{}:{}: Unknown command {name} at",start.row+1, start.column+1)]
    UnknownCommand { name: &'static str, start: Point },
    #[error("{}:{}: Unknown {kind}: {name} ",start.row+1, start.column+1)]
    UnknownCustom {
        kind: String,
        name: String,
        start: Point,
    },

    #[error("{0}")]
    ParseExpression(ParseExprError),
    #[error("{0}")]
    // TODO: This perhaps should not be a full error. Or make it incomplete node?
    PartialNode(String),
    #[error("syntax error")]
    InvalidSyntax,
}

#[derive(Debug, Error, Clone)]
#[error("An expected node was missing.")]
pub(crate) struct MissingNode;

impl From<MissingNode> for FromNodeError {
    fn from(_: MissingNode) -> Self {
        Self::PartialNode("Missing node".to_owned())
    }
}

impl<'a> IntoError<Node<'a>> for Option<Node<'a>> {
    type Error = MissingNode;
    fn into_err(self) -> Result<Node<'a>, Self::Error> {
        self.ok_or(MissingNode)
    }
}

impl diagnostics::Issue for SpannedError<FromNodeError> {
    fn diagnostic(&self) -> diagnostics::Diagnostic {
        let name = match self.error {
            FromNodeError::ParseIntError(_) => "parsing int failed",
            FromNodeError::ParseFloatError(_) => "parsing float failed",
            FromNodeError::UnknownCommand { .. } => "unknown command",
            FromNodeError::UnknownCustom { .. } => "custom error",
            FromNodeError::ParseExpression(_) => "expression error",
            FromNodeError::PartialNode(_) => "incomplete command",
            FromNodeError::InvalidSyntax => "invalid syntax",
        };

        diagnostics::Diagnostic {
            name: name.to_string(),
            severity: diagnostics::Severity::Error,
            span: self.span,
            message: self.to_string(),
        }
    }
}

impl From<ParseExprError> for FromNodeError {
    fn from(v: ParseExprError) -> Self {
        Self::ParseExpression(v)
    }
}

impl From<std::num::ParseFloatError> for FromNodeError {
    fn from(v: std::num::ParseFloatError) -> Self {
        Self::ParseFloatError(v)
    }
}

impl From<std::num::ParseIntError> for FromNodeError {
    fn from(v: std::num::ParseIntError) -> Self {
        Self::ParseIntError(v)
    }
}

// Traits that might make using `tree-sitter` more ergonomic. However, not sure.
// TODO: Remove these if they never get used.

#[allow(dead_code)]
trait GetChild {
    fn get_child(&self, index: usize) -> Result<Node, MissingChild>;
}

#[allow(dead_code)]
trait SeekChild {
    /// Moves the cursor to the location of `index` child and returns it as a result.
    fn seek_child(&mut self, index: usize) -> Result<Node, MissingChild>;
}

#[derive(Debug, Default, Error, Clone)]
#[error("Could not find child at index {0}")]
struct MissingChild(usize);

impl GetChild for Node<'_> {
    fn get_child(&self, index: usize) -> Result<Node, MissingChild> {
        self.child(index).ok_or(MissingChild(index))
    }
}

impl SeekChild for TreeCursor<'_> {
    fn seek_child(&mut self, index: usize) -> Result<Node, MissingChild> {
        if !self.goto_first_child() {
            return Err(MissingChild(index));
        }

        if index == 0 {
            return Ok(self.node());
        }

        for _ in 1..index {
            if !self.goto_next_sibling() {
                return Err(MissingChild(index));
            }
        }

        Ok(self.node())
    }
}

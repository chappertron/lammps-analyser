use tree_sitter::{Node, Point};

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
            name,
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

impl From<SpannedError<FromNodeError>> for FromNodeError {
    /// HACK: Allows ? on Spanned error for function that returns non-spanned.
    fn from(value: SpannedError<FromNodeError>) -> Self {
        value.error
    }
}

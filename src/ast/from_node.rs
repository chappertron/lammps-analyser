use owo_colors::OwoColorize;
use thiserror::Error;
use tree_sitter::{Node, Point, TreeCursor};

use crate::{diagnostic_report::ReportSimple, utils::into_error::IntoError};

use super::expressions::ParseExprError;

/// Converts from a tree-sitter node to the type this trait is implemented for.
///
/// Requires providing the text source, which is needed for extracting names and identifiers.
pub trait FromNode: Sized {
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, FromNodeError>;
}

// #[derive(Debug, Error, Clone)]
// #[error("Could not parse node from {start} to {end}: {kind}")]
// pub struct FromNodeError {
//     pub start: Point,
//     pub end: Point,
//     pub kind: FromNodeErrorType,
// }

// TODO: Just wrap expression parse error??
#[derive(Debug, Error, Clone)]
pub enum FromNodeError {
    #[error("Could not parse text as UTF-8 {0}")]
    Utf8Error(std::str::Utf8Error),
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
}

#[derive(Debug, Error, Clone)]
#[error("An expected node was missing.")]
pub(crate) struct MissingNode;

impl From<MissingNode> for FromNodeError {
    fn from(_: MissingNode) -> Self {
        Self::PartialNode("Missing node".to_string())
    }
}

impl<'a> IntoError<Node<'a>> for Option<Node<'a>> {
    type Error = MissingNode;
    fn into_err(self) -> Result<Node<'a>, Self::Error> {
        self.ok_or(MissingNode)
    }
}

impl ReportSimple for FromNodeError {
    fn make_simple_report(&self) -> String {
        match self {
            Self::Utf8Error(e) => format!("Could not parse text as UTF-8 {}", e.bright_red()),
            Self::ParseIntError(e) => format!("Could not parse text as int {}", e.bright_red()),
            Self::ParseFloatError(e) => format!("Could not parse text as float {}", e.bright_red()),
            Self::UnknownCommand { name, start } => format!(
                "{}:{}: Unknown command {} at",
                start.row + 1,
                start.column + 1,
                name.bright_red(),
            ),
            Self::UnknownCustom { kind, name, start } => format!(
                "{}:{}: Unknown {}: {} ",
                start.row + 1,
                start.column + 1,
                kind.bright_red(),
                name.bright_red(),
            ),

            Self::ParseExpression(e) => format!("{}", e.bright_red()),
            Self::PartialNode(e) => format!("{}", e.bright_red()),
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

impl From<std::str::Utf8Error> for FromNodeError {
    fn from(v: std::str::Utf8Error) -> Self {
        Self::Utf8Error(v)
    }
}

// Traits that might make using `tree-sitter` more ergonomic. However, not sure.

trait GetChild {
    fn get_child(&self, index: usize) -> Result<Node, MissingChild>;
}

trait SeekChild {
    /// Moves the cursor to the location of `index` child and resturns it as a result.
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

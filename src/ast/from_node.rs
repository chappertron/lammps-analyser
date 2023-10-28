use thiserror::Error;
use tree_sitter::{Node, TreeCursor};

trait FromNode: Sized {
    fn from_node(node: &Node, text: &[u8]) -> Result<Self, FromNodeError>;
}

// TODO: Just wrap expression parse error??
#[derive(Debug, Error, Clone)]
enum FromNodeError {
    #[error("Could not parse text as UTF-8 {0}")]
    Utf8Error(std::str::Utf8Error),
    #[error("Could not parse text as int {0}")]
    ParseIntError(std::num::ParseIntError),
    #[error("Could not parse text as float {0}")]
    ParseFloatError(std::num::ParseFloatError),
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

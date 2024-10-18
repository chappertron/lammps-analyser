use tree_sitter::Node;

use crate::ast::from_node::FromNodeError;

/// An extension trait that adds additional methods to `tree_sitter::Node`
pub(crate) trait NodeExt {
    /// A helper method for getting the text of a node from valid utf-8.
    ///
    /// # Panics
    /// Panics if node offset somehow is on invalid UTF-8 character boundaries.
    /// This should only happen if the file wasn't UTF-8 encoded.
    fn str_text<'a>(&self, source: &'a str) -> &'a str;

    /// The number of lines the node spans.
    fn n_lines(&self) -> usize;
}

impl NodeExt for Node<'_> {
    fn str_text<'a>(&self, source: &'a str) -> &'a str {
        // Should not panic as source is valid utf-8
        self.utf8_text(source.as_bytes()).expect("invalid utf-8")
    }

    fn n_lines(&self) -> usize {
        let tree_sitter::Range {
            start_point,
            end_point,
            ..
        } = self.range();

        end_point.row - start_point.row + 1
    }
}

pub trait ExpectNode {
    type Output;
    type Error;

    /// Convert an Optional<Node> or a Missing / Error node into an Error
    fn expect_node(self, message: impl Into<String>) -> Result<Self::Output, Self::Error>;
    fn expect_kind(
        self,
        kind: &str,
        message: impl Into<String>,
    ) -> Result<Self::Output, Self::Error>;
}

impl<'a> ExpectNode for Option<Node<'a>> {
    type Output = Node<'a>;
    type Error = FromNodeError;

    fn expect_node(self, message: impl Into<String>) -> Result<Node<'a>, Self::Error> {
        match self {
            Some(x) => {
                if x.is_error() | x.is_missing() {
                    return Err(FromNodeError::PartialNode(message.into()));
                }

                Ok(x)
            }
            None => Err(FromNodeError::PartialNode(message.into())),
        }
    }

    fn expect_kind(self, kind: &str, message: impl Into<String>) -> Result<Node<'a>, Self::Error> {
        match self {
            Some(x) => {
                if x.kind() != kind {
                    return Err(FromNodeError::PartialNode(message.into()));
                }

                Ok(x)
            }
            None => Err(FromNodeError::PartialNode(message.into())),
        }
    }
}

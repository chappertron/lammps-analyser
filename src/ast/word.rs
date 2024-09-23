use std::fmt::Display;

use tree_sitter::Node;

use crate::{spans::Span, utils::tree_sitter_helpers::NodeExt};

use super::from_node::{FromNode, FromNodeError};

#[derive(Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Word {
    pub contents: String,
    pub span: Span,
}
impl Display for Word {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.contents.fmt(f)
    }
}

impl Word {
    pub fn new(contents: String, span: impl Into<Span>) -> Self {
        let span = span.into();

        Self { contents, span }
    }

    pub fn as_str(&self) -> &str {
        self.contents.as_str()
    }

    pub(crate) fn parse_word(node: &Node, text: &str) -> Self {
        let contents = node.str_text(text.as_ref()).to_owned();
        let span = node.range().into();

        Word { contents, span }
    }
}

impl FromNode for Word {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: &str) -> Result<Self, Self::Error> {
        Ok(Word::parse_word(node, text))
    }
}

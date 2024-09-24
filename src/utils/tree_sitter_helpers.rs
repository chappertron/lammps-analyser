use tree_sitter::Node;

/// An extension trait that adds additional methods to `tree_sitter::Node`
pub(crate) trait NodeExt {
    /// A helper method for getting the text of a node from valid utf-8.
    ///
    /// # Panics
    /// Panics if node offset somehow is on invalid UTF-8 character boundaries.
    /// This should only happen if the file wasn't UTF-8 encoded.
    fn str_text<'a>(&self, source: &'a str) -> &'a str;
}

impl NodeExt for Node<'_> {
    fn str_text<'a>(&self, source: &'a str) -> &'a str {
        // Should not panic as source is valid utf-8
        self.utf8_text(source.as_bytes()).expect("invalid utf-8")
    }
}

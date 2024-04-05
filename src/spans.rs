/// A point within a file.
/// Zero-based, so the first row is zero.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Hash)]
pub struct Point {
    pub row: usize,
    pub column: usize,
}

impl From<tree_sitter::Point> for Point {
    fn from(value: tree_sitter::Point) -> Self {
        Point {
            row: value.row,
            column: value.column,
        }
    }
}

impl From<lsp_types::Position> for Point {
    fn from(value: lsp_types::Position) -> Self {
        Point {
            row: value.line as usize,
            column: value.character as usize,
        }
    }
}

impl Point {
    /// Converts into a `tree_sitter::Point`
    pub fn into_tree_sitter(self) -> tree_sitter::Point {
        tree_sitter::Point {
            row: self.row,
            column: self.column,
        }
    }

    /// Converts into an `lsp_types::Position`
    /// Note: may cause issues if row/column > `u32::MAX`
    pub fn into_lsp_type(self) -> lsp_types::Position {
        lsp_types::Position {
            line: self.row as u32,
            character: self.column as u32,
        }
    }
}

#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Hash)]
pub struct Span {
    pub start: Point,
    pub end: Point,
}

impl From<tree_sitter::Range> for Span {
    fn from(value: tree_sitter::Range) -> Self {
        Self {
            start: value.start_point.into(),
            end: value.end_point.into(),
        }
    }
}

impl From<lsp_types::Range> for Span {
    fn from(value: lsp_types::Range) -> Self {
        Self {
            start: value.start.into(),
            end: value.end.into(),
        }
    }
}

impl Span {
    /// Convert the span into a `tree-sitter::Range`.
    /// Requires additional information about the bytes span of the token to be provided.
    pub fn into_tree_sitter(self, start_byte: usize, end_byte: usize) -> tree_sitter::Range {
        tree_sitter::Range {
            start_point: self.start.into_tree_sitter(),
            end_point: self.end.into_tree_sitter(),
            start_byte,
            end_byte,
        }
    }

    /// Convert the span into an `lsp_types::Range`
    /// NOTE: if the line/column number is > `u32::MAX`, casting from `usize` may result in issues.
    pub fn into_lsp_types(self) -> lsp_types::Range {
        lsp_types::Range {
            start: self.start.into_lsp_type(),
            end: self.end.into_lsp_type(),
        }
    }
}

impl Span {
    pub fn in_range(&self, range: &Span) -> bool {
        range.start <= self.start && self.end <= range.end
    }
}

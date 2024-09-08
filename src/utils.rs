use crate::spans::Point;
use lsp_types::Position;

use crate::identifinder::{IdentiFinder, NameAndType};

/// Convert a tree-sitter point to an lsp-position
pub fn point_to_position(point: &tree_sitter::Point) -> lsp_types::Position {
    Position {
        line: point.row as u32,
        character: point.column as u32,
    }
}

/// Convert a lsp-position to an tree-sitter point
pub fn position_to_point(pos: &lsp_types::Position) -> tree_sitter::Point {
    tree_sitter::Point {
        row: pos.line as usize,
        column: pos.character as usize,
    }
}

/// convert a `tree_sitter::Range` into a `lsp_types::Range`
/// Note, the reverse is not trivial because a `tree_sitter::Range` requires byte information
/// that is not used in `lsp_types`
pub fn ts_range_to_lsp_range(range: &tree_sitter::Range) -> lsp_types::Range {
    lsp_types::Range {
        start: point_to_position(&range.start_point),
        end: point_to_position(&range.end_point),
    }
}

/// Find the symbol under the given point in the file
pub fn get_symbol_at_point<'a>(
    identifinder: &'a IdentiFinder,
    point: &Point,
) -> Option<&'a NameAndType> {
    let symbols = identifinder.symbols();

    // Might've created something a little slow here...
    for (k, v) in symbols {
        for r in v.refs().iter() {
            if *point >= r.start() && *point <= r.end() {
                // TODO: Return multiple if they exist here????
                return Some(k);
            }
        }
    }

    None
}

pub trait InTSRange {
    /// Returns true if point is entirely in the `Range`
    fn in_range(&self, range: &tree_sitter::Range) -> bool;
}

impl InTSRange for tree_sitter::Point {
    fn in_range(&self, range: &tree_sitter::Range) -> bool {
        &range.start_point <= self && self <= &range.end_point
    }
}

impl InTSRange for tree_sitter::Range {
    fn in_range(&self, range: &tree_sitter::Range) -> bool {
        range.start_point <= self.start_point && self.end_point <= range.end_point
    }
}

/// TODO: Perhaps mixing TS and LSP points/positions will get confusing...
impl InTSRange for lsp_types::Position {
    fn in_range(&self, range: &tree_sitter::Range) -> bool {
        let ts_point = position_to_point(self);
        range.start_point <= ts_point && ts_point <= range.end_point
    }
}

pub(crate) mod into_error {

    /// A trait designed for allowing easier conversions from `Option<T>` to `Result<T, E>`.
    ///
    /// This needed for foreign types T to circumvent the orphan rules, that prevent using the `From` trait.
    /// This can only be defined Once per type, but then From<Error> can be implemented for other
    /// error types with the question mark operator.
    pub(crate) trait IntoError<T> {
        type Error;
        // TODO: allow for custom messages
        fn into_err(self) -> Result<T, Self::Error>;
    }
}

pub(crate) mod parsing {
    use tree_sitter::Parser;

    pub fn setup_parser() -> Parser {
        let mut parser = Parser::new();
        parser.set_language(tree_sitter_lammps::language()).unwrap();
        parser
    }
}

#[cfg(test)]
pub(crate) mod testing {
    use tree_sitter::Parser;

    pub fn setup_parser() -> Parser {
        super::parsing::setup_parser()
    }
    pub(crate) fn parse(source_bytes: impl AsRef<[u8]>) -> tree_sitter::Tree {
        let source_bytes = source_bytes.as_ref();
        let mut parser = setup_parser();
        parser.parse(source_bytes, None).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn in_range_equal() {
        use tree_sitter::{Point, Range};
        let row_0 = 1;
        let row_1 = 1;
        let column_0 = 0;
        let column_1 = 20;

        let range_a = Range {
            start_point: Point {
                row: row_0,
                column: column_0,
            },
            start_byte: 0,
            end_byte: 0,
            end_point: Point {
                row: row_1,
                column: column_1,
            },
        };

        let range_b = Range {
            start_point: Point {
                row: row_0,
                column: column_0,
            },
            start_byte: 10, // Different offset, because shouldn't affect range check
            end_byte: 10,
            end_point: Point {
                row: row_1,
                column: column_1,
            },
        };

        assert!(range_a.in_range(&range_b));
        assert!(range_b.in_range(&range_a));
    }
}

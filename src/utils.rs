use lsp_types::Position;
use tree_sitter::Point;

use crate::identifinder::{self, Ident, IdentiFinder, NameAndType};

/// Convert a tree-sitter point to an lsp-position
pub fn point_to_position(point: &Point) -> Position {
    Position {
        line: point.row as u32,
        character: point.column as u32,
    }
}

/// Find the symbol at the given point in the file
/// TODO Return multiple if they exist here????
pub fn get_symbol_at_point(identifinder: &IdentiFinder, point: Point) -> Option<&NameAndType> {
    let symbols = identifinder.symbols();

    // Might've created something a little slow here...
    for (k, v) in symbols {
        for r in v.refs().iter() {
            if point >= r.start && point <= r.end {
                return Some(k);
            }
        }
    }

    None
}

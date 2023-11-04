use lsp_types::Position;
use tree_sitter::Point;

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
    Point {
        row: pos.line as usize,
        column: pos.character as usize,
    }
}

/// Find the symbol at the given point in the file
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

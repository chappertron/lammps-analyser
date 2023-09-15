use lsp_types::Position;
use tree_sitter::Point;

/// Convert a tree-sitter point to an lsp-position
pub fn point_to_position(point: &Point) -> Position {
    Position {
        line: point.row as u32,
        character: point.column as u32,
    }
}

use anyhow::{Ok, Result};
use owo_colors::OwoColorize;
use thiserror::Error;
use tree_sitter::{ Point, Query, QueryCursor, Tree, TreeCursor};
pub struct ErrorFinder {
    pub query: Query,
    cursor: QueryCursor,
    syntax_errors: Vec<SyntaxError>,
}

impl ErrorFinder {
    pub fn new() -> Result<Self> {
        let query = Query::new(
            tree_sitter_lammps::language(),
            "
            (ERROR) @syntax_error
            ;(MISSING) @syntax_error
            ",
        )?;
        let cursor = QueryCursor::new();
        let syntax_errors = vec![];

        Ok(Self {
            query,
            cursor,
            syntax_errors,
        })
    }

    pub fn syntax_errors(&self) -> &[SyntaxError] {
        self.syntax_errors.as_ref()
    }

    pub fn find_syntax_errors(
        &mut self,
        tree: &Tree,
        source_code: &[u8],
    ) -> Result<&[SyntaxError]> {
        let matches = self
            .cursor
            .matches(&self.query, tree.root_node(), source_code);
        for mat in matches {
            let text = mat.captures[0].node.utf8_text(source_code)?;
            let start = mat.captures[0].node.start_position();
            let end = mat.captures[0].node.end_position();
            self.syntax_errors.push(SyntaxError {
                text: text.into(),
                start,
                end,
            });
        }
        Ok(self.syntax_errors())
    }

    /// Tree-sitter can't currently query for missing nodes, so recursivley walking the tree instead
    /// Missing nodes are also not reported as errors, so this is needed?
    /// TODO Walk back from missing nodes to work out the proper node? 
    pub fn find_missing_nodes(
        &mut self,
        tree: &Tree,
    ) -> Result<&Vec<SyntaxError>> {
        let mut cursor = tree.root_node().walk();
        let mut missing_nodes = vec![];
        fn recur_missing(
            cursor: &mut TreeCursor,
            missing_nodes: &mut Vec<(Point, Point)>,
        ) {
            if cursor.node().child_count() == 0 {
                if cursor.node().is_missing() {
                    // println!("{} {}:{}","Missing Node:".red(),cursor.node().start_position().row+1,cursor.node().start_position().column+1);
                    let node = cursor.node();
                    missing_nodes.push((
                        node.start_position(),
                        node.end_position(),
                    ));
                }
            } else {
                // Go to the first child, then recur
                cursor.goto_first_child();
                recur_missing(cursor, missing_nodes);
                // Go to the  siblings
                while cursor.goto_next_sibling() {
                    recur_missing(cursor, missing_nodes);
                }
                cursor.goto_parent();
            }
        }

        recur_missing(&mut cursor, &mut missing_nodes);
        self.syntax_errors.extend(missing_nodes
            .iter()
            .map(|(start_position, end_position)| {
                SyntaxError {
                    text: "Missing Node".to_string(),
                    start: *start_position,
                    end: *end_position,
                }
            })
            );

        Ok(&self.syntax_errors)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Hash, Error)]
#[error("{}:{}: {} `{}`",start.row+1,start.column+1,"Invalid Syntax".bright_red(),text)]
pub struct SyntaxError {
    pub text: String,
    pub start: Point,
    pub end: Point,
}

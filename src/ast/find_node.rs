use crate::spans::Point;

use super::{Ast, Command};

impl Ast {
    /// Find the Command node located at a given text location.
    ///
    /// Uses binary search over the nodes as these should be in order of span.
    #[allow(clippy::missing_panics_doc)] // expect shouldn't panic here.
    pub fn find_node(&self, point: Point) -> Option<&Command> {
        self.commands
            .binary_search_by(|node| {
                // NOTE: This expect is OK. partial cmp between `Point` and `Span` never returns `None`
                #[allow(clippy::expect_used)]
                node.span()
                    .partial_cmp(&point)
                    // Invert. this comparison says whether the point is in the range or not.
                    .expect("This should always return `Some`")
            })
            .map(|i| dbg!(i))
            .map_err(|e| dbg!(e))
            .ok()
            .map(|i| &self.commands[i])
    }
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {

    use crate::{
        ast::{Command, FixDef},
        styles::FixStyle,
    };

    use super::*;

    #[test]
    fn find_node_test() {
        let text = r#"# this is a comment


fix     myfix all nve
compute mycompute all temp 
compute mycompute2 all ke/atom
"#;

        dbg!(&text);

        let n_lines = text.lines().count();

        let mut parser = crate::utils::testing::setup_parser();
        let tree = parser.parse(text, None).expect("Should be able to parse");

        let ast = super::super::ts_to_ast(&tree, text).expect("Failed to generated AST");
        dbg!(&ast.commands);

        // Don't care about exact definition, just the type
        let is_fix_def = match ast.find_node(Point { row: 3, column: 0 }) {
            Some(
                Command::Fix(FixDef {
                    fix_style: FixStyle::Nve,
                    ..
                }),
                ..,
            ) => true,
            node => {
                dbg!(node);
                false
            }
        };

        assert!(is_fix_def);

        // Out of range
        assert_eq!(
            ast.find_node(Point {
                row: n_lines + 10,
                column: 0
            }),
            None
        );
    }
}

//! Convert the treesitter trees into an AST.

// For denying unwraps and expects in this file
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]

// Sub modules.
pub mod find_node;
pub mod from_node;
pub mod impl_helpers;

mod arguments;
mod commands;
mod compute_def;
mod expressions;
mod fix_def;
mod ident;
mod variable_def;
mod word;

// Re-exports
pub use arguments::{Argument, ArgumentKind};
pub use commands::{Command, GenericCommand};
pub use compute_def::ComputeDef;
pub use expressions::Expression;
pub use fix_def::FixDef;
pub use variable_def::VariableDef;
pub use word::Word;

pub use ident::{Ident, IdentType};

use crate::spanned_error::SpannedError;

use self::from_node::{FromNode, FromNodeError};
use tree_sitter::Tree;

#[derive(Debug, Clone, PartialEq, Default)]
/// A list of commands
pub struct Ast {
    pub(crate) commands: Vec<Command>,
}

/// An AST that could not be fully parsed. Also contains errors.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialAst {
    pub ast: Ast,
    pub errors: Vec<SpannedError<FromNodeError>>,
}

/// Transform a `tree-sitter::Tree` into an Abstract Syntax Tree
pub fn ts_to_ast(tree: &Tree, text: &str) -> Result<Ast, PartialAst> {
    let mut cursor = tree.walk();

    let mut commands = vec![];
    let mut errors = vec![];
    cursor.goto_first_child();

    for node in tree.root_node().named_children(&mut cursor) {
        // TODO: This if-statement may be removable
        if node.kind() != "comment" {
            let result = Command::from_node(&node, text); //.with_span(node.range().into());

            match result {
                Ok(cmd) => {
                    commands.push(cmd);
                }
                Err((cmd, error)) => {
                    // Add ERROR nodes.
                    commands.push(cmd);
                    errors.push(SpannedError::new(error, node.range()));
                }
            }
        }
    }

    if errors.is_empty() {
        Ok(Ast { commands })
    } else {
        Err(PartialAst {
            ast: Ast { commands },
            errors,
        })
    }
}

#[cfg(test)]
// Allow unwraps in the tests module, but not in the parent module.
#[allow(clippy::unwrap_used)]
#[allow(clippy::expect_used)]
mod tests {

    use serde_json::Value;

    // use pretty_assertions::assert_eq;
    use crate::utils::testing::parse;

    use super::ts_to_ast;

    #[test]
    #[ignore = "Not yet fully implemented"]
    fn test_ast() {
        // let source_bytes = include_bytes!("../../fix.lmp");
        let source = include_str!("../../example_input_scripts/in.nemd");
        let tree = parse(source);

        let _ast = ts_to_ast(&tree, source);
        // dbg!(ast.unwrap());

        unimplemented!()
    }

    #[test]
    fn node_coverage() {
        const NODES_STR: &str = tree_sitter_lammps::NODE_TYPES;

        std::fs::write("./node_types.json", NODES_STR).expect("Failed to write node types out.");
        let value: Value = serde_json::from_str(NODES_STR).expect("Failed to deserialise");

        dbg!(&value);

        // dbg!(NODES_STR);
        for node_type in value.as_array().expect("top level should be an array") {
            dbg!(node_type);
        }

        unimplemented!()
    }
}

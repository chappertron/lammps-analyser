use tree_sitter::Node;

use crate::spans::{Point, Span};

use super::from_node::{FromNode, FromNodeError};
use crate::ast;

#[derive(Debug, PartialEq, Clone)]
pub enum Command {
    Generic(GenericCommand),
    /// A Fix definition
    Fix(ast::FixDef),
    /// A compute definition
    Compute(ast::ComputeDef),
    /// A variable definition
    VariableDef(ast::VariableDef),
    Shell(Span),
    Error(Span),
}

impl Command {
    pub fn span(&self) -> Span {
        match self {
            Command::Generic(cmd) => cmd.span(),
            Command::Fix(cmd) => cmd.span,
            Command::Compute(cmd) => cmd.span,
            Command::VariableDef(cmd) => cmd.span,
            Command::Shell(span) => *span,
            Command::Error(span) => *span,
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct GenericCommand {
    pub name: ast::Word,
    pub args: Vec<ast::Argument>,
    pub start: Point,
    pub end: Point,
    pub start_byte: usize,
    pub end_byte: usize,
}

impl GenericCommand {
    pub fn span(&self) -> Span {
        Span {
            start: self.start,
            end: self.end,
        }
    }
}

impl FromNode for GenericCommand {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: &str) -> Result<Self, FromNodeError> {
        let mut cursor = node.walk();
        let start = node.start_position();
        let end = node.end_position();
        let start_byte = node.start_byte();
        let end_byte = node.end_byte();

        let mut args = vec![];

        debug_assert!(cursor.node() == *node);

        if !cursor.goto_first_child() {
            return Err(FromNodeError::PartialNode(
                "missing command name".to_owned(),
            ));
        }

        // TODO: use a field in the TS grammar instead?
        let name = ast::Word::from_node(&cursor.node(), text)?;

        while cursor.goto_next_sibling() {
            for node in cursor.node().children(&mut cursor) {
                args.push(ast::Argument::from_node(&node, text)?);
            }
        }

        Ok(GenericCommand {
            name,
            args,
            start: start.into(),
            end: end.into(),
            start_byte,
            end_byte,
        })
    }
}

impl FromNode for Command {
    // NOTE: This is chosen as the error type so an error node can be returned instead
    type Error = (Self, FromNodeError);
    fn from_node(node: &Node, text: &str) -> Result<Self, (Self, FromNodeError)> {
        let error_span = |e| (Self::Error(node.range().into()), e);
        let mut result = match node.kind() {
            "fix" => Ok(Self::Fix(
                ast::FixDef::from_node(node, text).map_err(error_span)?,
            )),
            "compute" => Ok(Self::Compute(
                ast::ComputeDef::from_node(node, text).map_err(error_span)?,
            )),
            // TODO: Make a variable deltion it's own type
            "variable_def" | "variable_del" => Ok(Self::VariableDef(
                ast::VariableDef::from_node(node, text).map_err(error_span)?,
            )),
            "shell" => Ok(Self::Shell(node.range().into())),
            // Fall back to the generic command type
            "command" => Ok(Self::Generic(
                GenericCommand::from_node(node, text).map_err(error_span)?,
            )),
            "ERROR" => Ok(Self::Error(node.range().into())),

            // NOTE: make this variant a panic for testing purposes.
            #[cfg(feature = "ast_panics")]
            c => panic!("unknown command kind {c}"),
            _ => Ok(Self::Error(node.range().into())),
        };

        if let Ok(Self::Error(span)) = result {
            // NOTE: Check the child node and infer what type of node it was supposed to be
            // and pass to that parser.

            result = match node.child(0).map(|node| node.kind()) {
                // Try and parse this as a compute
                Some("compute") => Ok(Self::Compute(
                    ast::ComputeDef::from_node(node, text).map_err(error_span)?,
                )),
                Some("fix_id") => Ok(Self::Fix(
                    ast::FixDef::from_node(node, text).map_err(error_span)?,
                )),
                Some("variable") => Ok(Self::VariableDef(
                    ast::VariableDef::from_node(node, text).map_err(error_span)?,
                )),
                // Cannot further process
                // NOTE: This is `Ok` rather than `Error` because syntax errors are also
                // found within `ErrorFinder`
                _ => Ok(Self::Error(span)),
            };
        }

        result
    }
}

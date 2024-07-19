//! Convert the treesitter trees into an AST.

// For denying unwraps and expects in this file
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]

pub mod expressions;
pub mod find_node;
pub mod from_node;

use crate::{
    commands::CommandName,
    spanned_error::{SpannedError, WithSpan},
    spans::{Point, Span},
};
use std::fmt::Display;

use tree_sitter::{Node, Tree};

use crate::{
    compute_styles::ComputeStyle, fix_styles::FixStyle, identifinder::Ident,
    utils::into_error::IntoError,
};

use self::from_node::{FromNode, FromNodeError};

#[derive(Debug, Clone, PartialEq, Default)]
/// A list of commands
/// Perhaps not truly an AST
pub struct Ast {
    pub commands: Vec<CommandNode>,
}

/// An AST that could not be fully produced. Also contains errors.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialAst {
    pub ast: Ast,
    pub errors: Vec<SpannedError<FromNodeError>>,
}

impl Ast {
    /// Find the command corresponding to a given point
    pub fn find_point(&self, point: &Point) -> Option<&CommandNode> {
        self.commands
            .iter()
            .find(|cmd| cmd.range.start <= *point && cmd.range.end >= *point)
    }
}

/// TODO: return both the partial AST and the Errors
pub fn ts_to_ast(tree: &Tree, text: impl AsRef<[u8]>) -> Result<Ast, PartialAst> {
    let mut cursor = tree.walk();
    let text = text.as_ref();

    let mut commands = vec![];
    let mut errors = vec![];
    cursor.goto_first_child();

    loop {
        // Advance cursor and skip if a comment
        if cursor.goto_first_child() && cursor.node().kind() != "comment" {
            let result = CommandNode::from_node(&cursor.node(), text)
                .with_span(cursor.node().range().into());

            match result {
                Ok(cmd) => commands.push(cmd),
                Err(err) => errors.push(err),
            }

            cursor.goto_parent();
        }

        // If no more commands, break!
        if !cursor.goto_next_sibling() {
            break;
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

/// A command in the LAMMPS Input Script
///  TODO: rename this to be a command node???
#[derive(Debug, PartialEq, Clone)]
pub struct CommandNode {
    pub command_type: CommandType,
    range: Span,
}

impl CommandNode {
    pub fn range(&self) -> Span {
        self.range
    }
}

impl FromNode for CommandNode {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        let range = node.range().into();

        let command_type = CommandType::from_node(node, text)?;

        Ok(CommandNode {
            command_type,
            range,
        })
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum CommandType {
    GenericCommand(GenericCommand),
    NamedCommand(NamedCommand),
}

impl FromNode for CommandType {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, FromNodeError> {
        let cmd = if NamedCommand::try_from(node.kind()).is_ok() {
            // TODO: add arguments
            CommandType::NamedCommand(NamedCommand::from_node(node, text)?)
        } else {
            CommandType::GenericCommand(GenericCommand::from_node(node, text)?)
        };

        Ok(cmd)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct GenericCommand {
    pub name: CommandName,
    pub args: Vec<Argument>,
    pub start: Point,
    pub end: Point,
    pub start_byte: usize,
    pub end_byte: usize,
}

impl FromNode for GenericCommand {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, FromNodeError> {
        let mut cursor = node.walk();
        // let kind = node.kind().to_string();
        let start = node.start_position();
        let end = node.end_position();
        let start_byte = node.start_byte();
        let end_byte = node.end_byte();

        let text = text.as_ref();

        let mut args = vec![];

        debug_assert!(cursor.node() == *node);

        cursor.goto_first_child();

        // TODO: use a field in the TS grammar
        let name = cursor.node().utf8_text(text)?;

        while cursor.goto_next_sibling() {
            for node in cursor.node().children(&mut cursor) {
                args.push(Argument::from_node(&node, text)?);
            }
        }

        cursor.goto_parent();
        Ok(GenericCommand {
            name: name.into(),
            args,
            start: start.into(),
            end: end.into(),
            start_byte,
            end_byte,
        })
    }
}

/// Acceptable argument types for LAMMPS commands
#[derive(Debug, PartialEq, Clone)]
pub enum Argument {
    /// Expression from LAMMPS
    Int(isize),
    Float(f64),
    Bool(bool),
    ArgName(String), // TODO: Rename to keyword arg?
    /// Variables within curly braces
    VarCurly(Ident),
    VarRound(expressions::Expression),
    String(String),
    Expression(expressions::Expression),
    // TODO: Remove? Can't know if a group name until further on in the process???
    // Perhaps make it an identifier that then is decided to be either
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    UnderscoreIdent(Ident),
    /// A white-space delimited word
    Word(String),
    /// Indicates an invalid node in the tree-sitter grammar
    Error,
}

impl FromNode for Argument {
    type Error = FromNodeError;
    fn from_node(
        node: &Node,
        // _cursor: &mut tree_sitter::TreeCursor,
        text: impl AsRef<[u8]>,
    ) -> Result<Self, <Argument as FromNode>::Error> {
        // TODO: make these variants more complete.
        //.child(0).unwrap()
        // Did removing above from match fix things???
        let text = text.as_ref();
        match node.kind() {
            "int" => Ok(Self::Int(node.utf8_text(text)?.parse::<isize>()?)),
            "float" => Ok(Self::Float(node.utf8_text(text)?.parse::<f64>()?)),
            "bool" => Ok(Self::Bool(match node.utf8_text(text)? {
                "on" | "yes" | "true" => true,
                "off" | "no" | "false" => false,
                invalid_name => Err(FromNodeError::UnknownCustom {
                    kind: "Bool Variant".to_owned(),
                    name: invalid_name.to_string(),
                    start: node.start_position(),
                })?,
            })),
            // TODO: Expressions not wrapped in varround are not valid???
            "expression" => Ok(Self::Expression(expressions::Expression::parse_expression(
                node, text,
            )?)),
            "string" => Ok(Self::String(node.utf8_text(text)?.to_owned())),
            "group" => Ok(Self::Group),
            "underscore_ident" => Ok(Self::UnderscoreIdent(Ident::new(
                &node.child(0).into_err()?,
                text,
            )?)),
            "var_curly" => Ok(Self::VarCurly(Ident::new(
                &node.child(1).into_err()?,
                text,
            )?)),
            "var_round" => Ok(Self::VarRound(expressions::Expression::parse_expression(
                &node.child(1).into_err()?,
                text,
            )?)),
            "argname" => Ok(Self::ArgName(
                //.child(0).into_err()?
                node.utf8_text(text)?.to_string(),
            )),
            "word" => Ok(Self::Word(
                //.child(0).into_err()?
                node.utf8_text(text)?.to_string(),
            )),
            // TODO: It seems weird sending something called error through ok.
            "ERROR" => Ok(Self::Error),
            x => Err(FromNodeError::UnknownCustom {
                kind: "argument type".to_string(),
                name: x.to_string(),
                start: node.start_position(),
            }),
        }
    }
}

impl Display for Argument {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Argument::Int(x) => write!(f, "int: {x}"),
            Argument::Float(x) => write!(f, "float: {x}"),
            Argument::Bool(x) => write!(f, "bool: {x}"),
            Argument::ArgName(x) => write!(f, "argname: {x}"),
            Argument::VarCurly(x) => write!(f, "var_curly: ${{{0}}}", x.name),
            Argument::VarRound(x) => write!(f, "var_round: $({x})"),
            // TODO: properly implement string
            Argument::String(s) => write!(f, "string: {s}"),
            Argument::Expression(x) => write!(f, "expression: {x}"),
            Argument::Group => write!(f, "group"),
            Argument::UnderscoreIdent(x) => write!(f, "underscore_ident: {x}"),
            Argument::Word(x) => write!(f, "word: {x}"),
            Argument::Error => write!(f, "ERROR"),
        }
    }
}

/// Commands that have a special form in the tree sitter grammar
/// TODO: add arguments
/// TODO: Add command location
#[derive(Debug, PartialEq, Clone)]
pub enum NamedCommand {
    Fix(FixDef),
    Compute(ComputeDef),
    Style,
    Modify,
    AtomStyle,
    Boundary,
    VariableDef,
    ThermoStyle,
    Thermo,
    Units,
    Run,
    Shell,
}

impl FromNode for NamedCommand {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<NamedCommand, Self::Error> {
        match node.kind() {
            "fix" => Ok(NamedCommand::Fix(FixDef::from_node(node, text)?)),
            "compute" => Ok(NamedCommand::Compute(ComputeDef::from_node(node, text)?)),
            "style" => Ok(NamedCommand::Style),
            "modify" => Ok(NamedCommand::Modify),
            "atom_style" => Ok(NamedCommand::AtomStyle),
            "boundary" => Ok(NamedCommand::Boundary),
            "variable" => Ok(NamedCommand::VariableDef),
            "thermo_style" => Ok(NamedCommand::ThermoStyle),
            "thermo" => Ok(NamedCommand::Thermo),
            "units" => Ok(NamedCommand::Units),
            "run" => Ok(NamedCommand::Run),
            "shell" => Ok(NamedCommand::Shell),
            _ => Err(FromNodeError::UnknownCommand {
                name: node.kind(),
                start: node.start_position(),
            }),
        }
    }
}

#[derive(Debug, Default, PartialEq, Clone)]
pub struct FixDef {
    pub fix_id: Ident,    // Or just keep as a string?
    pub group_id: String, // TODO:  Create group identifiers
    pub fix_style: FixStyle,
    // Arguments for fix command
    pub args: Vec<Argument>,
}

impl FixDef {
    pub fn range(&self) -> Span {
        self.fix_id.range()
    }
}

impl FromNode for FixDef {
    type Error = FromNodeError;

    /// TODO: Hand a cursor instead???
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        let mut cursor = node.walk();
        let text = text.as_ref();

        // TODO: handle case this is false
        cursor.goto_first_child();
        let mut children = node.children(&mut cursor);

        // skip the fix keyword
        children.next();
        let fix_id = Ident::new(
            &children
                .next()
                // TODO: This causes issues further up the chain.
                .ok_or(FromNodeError::PartialNode("Missing identifier".into()))?,
            text,
        )?;

        let group_id = children.next().into_err()?.utf8_text(text)?.into();
        let fix_style = children.next().into_err()?.utf8_text(text)?.into();

        let args = if let Some(args) = children.next() {
            // No longer needed beyond args. Lets us use cursor again
            drop(children);
            args.children(&mut cursor)
                .map(|x| Argument::from_node(&x, text))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            vec![]
        };

        Ok(FixDef {
            fix_id,
            group_id,
            fix_style,
            args,
        })
    }
}

#[derive(Debug, Default, PartialEq, Clone)]
pub struct ComputeDef {
    pub compute_id: Ident, // Or just keep as a string?
    pub group_id: String,  // TODO:  Create group identifiers
    pub compute_style: ComputeStyle,
    // Arguments for fix command
    pub args: Vec<Argument>,
}

impl ComputeDef {
    pub fn range(&self) -> Span {
        self.compute_id.range()
    }
}

impl FromNode for ComputeDef {
    type Error = FromNodeError;
    /// TODO: Hand a cursor instead???
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        let mut cursor = node.walk();

        // TODO: handle case this is false
        cursor.goto_first_child();
        let mut children = node.children(&mut cursor);
        let text = text.as_ref();

        // skip the fix keyword
        children.next();
        let compute_id = Ident::new(&children.next().into_err()?, text)?;
        let group_id = children.next().into_err()?.utf8_text(text)?.into();
        let compute_style = children.next().into_err()?.utf8_text(text)?.into();

        let args = if let Some(args) = children.next() {
            // No longer needed beyond args. Lets us use cursor again
            drop(children);
            args.children(&mut cursor)
                .map(|x| Argument::from_node(&x, text))
                .collect::<Result<Vec<_>, _>>()?
        } else {
            vec![]
        };

        Ok(ComputeDef {
            compute_id,
            group_id,
            compute_style,
            args,
        })
    }
}

/// This doesn't work for the new case that the fix has data
impl TryFrom<&str> for NamedCommand {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "fix" => Ok(Self::Fix(FixDef::default())),
            "compute" => Ok(Self::Compute(ComputeDef::default())),
            "style" => Ok(Self::Style),
            "modify" => Ok(Self::Modify),
            "atom_style" => Ok(Self::AtomStyle),
            "boundary" => Ok(Self::Boundary),
            "variable" => Ok(Self::VariableDef),
            "thermo_style" => Ok(Self::ThermoStyle),
            "thermo" => Ok(Self::Thermo),
            "units" => Ok(Self::Units),
            "run" => Ok(Self::Run),
            s => Err(format!("Unknown named command: {s}")),
        }
    }
}

#[cfg(test)]
// Allow unwraps in the tests module, but not in the parent module.
#[allow(clippy::unwrap_used)]
#[allow(clippy::expect_used)]
mod tests {
    use tree_sitter::Parser;

    use crate::ast::expressions::{BinaryOp, Expression};
    use crate::ast::from_node::FromNode;
    use crate::ast::{Argument, ComputeDef, FixDef};
    use crate::compute_styles::ComputeStyle;
    use crate::fix_styles::FixStyle;
    use crate::identifinder::{Ident, IdentType};

    use super::ts_to_ast;

    fn setup_parser() -> Parser {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        parser
    }

    #[test]
    #[ignore = "Not yet fully implemented"]
    fn test_ast() {
        let mut parser = setup_parser();
        // let source_bytes = include_bytes!("../../fix.lmp");
        let source_bytes = include_bytes!("../../example_input_scripts/in.nemd");
        let tree = parser.parse(source_bytes, None).unwrap();

        let ast = ts_to_ast(&tree, source_bytes);
        // dbg!(ast.unwrap());

        unimplemented!()
    }

    #[test]
    fn parse_fix_no_args() {
        let mut parser = setup_parser();
        let source_bytes = b"fix NVE all nve";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();
        dbg!(command_node.to_sexp());
        // Lots of tedium to parsing this...
        let fix_node = dbg!(command_node.child(0)).unwrap();

        dbg!(fix_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            FixDef::from_node(&fix_node, source_bytes.as_slice()).unwrap(),
            FixDef {
                fix_id: Ident {
                    name: "NVE".into(),
                    ident_type: IdentType::Fix,
                    ..Default::default()
                },
                group_id: "all".into(),
                fix_style: FixStyle::Nve,
                args: vec![],
            }
        );
    }

    #[test]
    #[ignore = "Not yet fully implemented"]
    // TODO: Finish this test
    fn parse_index_variable() {
        let mut parser = setup_parser();
        let text = b"variable index step4.1.aatm";

        let tree = parser.parse(text, None).unwrap();

        let ast = ts_to_ast(&tree, text);

        if ast.is_err() {
            dbg!(ast.unwrap_err());
        } else {
            dbg!(ast.unwrap());
        }

        unimplemented!()
    }

    #[test]
    fn parse_fix_with_args() {
        let mut parser = setup_parser();
        let source_bytes = b"fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();
        dbg!(command_node.to_sexp());
        // Lots of tedium to parsing this...
        let fix_node = dbg!(command_node.child(0)).unwrap();

        dbg!(fix_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            FixDef::from_node(&fix_node, source_bytes.as_slice()).unwrap(),
            FixDef {
                fix_id: Ident {
                    name: "NVT".into(),
                    ident_type: IdentType::Fix,
                    ..Default::default()
                },
                group_id: "all".into(),
                fix_style: FixStyle::Nvt,
                args: vec![
                    Argument::Word("temp".into()),
                    Argument::Int(1),
                    Argument::Float(1.5),
                    Argument::VarRound(Expression::BinaryOp(
                        Expression::Float(100.0).into(),
                        BinaryOp::Multiply,
                        Expression::ThermoKeyword("dt".into()).into(),
                    ))
                ],
            }
        );
    }

    #[test]
    fn parse_compute_with_args() {
        let mut parser = setup_parser();
        let source_bytes = b"compute T_hot water temp/region hot_region";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();
        dbg!(command_node.to_sexp());
        // Lots of tedium to parsing this...
        let compute_node = dbg!(command_node.child(0)).unwrap();

        dbg!(compute_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            ComputeDef::from_node(&compute_node, source_bytes.as_slice()).unwrap(),
            ComputeDef {
                compute_id: Ident {
                    name: "T_hot".into(),
                    ident_type: IdentType::Compute,
                    ..Default::default()
                },
                group_id: "water".into(),
                compute_style: ComputeStyle::TempRegion,
                // TODO: Change to a more generic word argument
                args: vec![Argument::Word("hot_region".into()),],
            }
        );
    }
}

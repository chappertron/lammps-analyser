//! Convert the treesiterr trees into an AST.

pub mod expressions;
pub mod from_node;
use std::fmt::Display;

use tree_sitter::{Node, Point, Range, Tree};

use crate::{fix_styles::FixStyle, identifinder::Ident};

use self::from_node::{FromNode, FromNodeError};

#[derive(Debug, Clone, PartialEq)]
/// A list of commands
/// Perhaps not truly an AST
pub struct Ast {
    pub commands: Vec<CommandNode>,
}

pub fn ts_to_ast(tree: &Tree, text: &[u8]) -> Result<Ast, FromNodeError> {
    let mut cursor = tree.walk();

    let mut commands = vec![];
    // println!("{}", tree.root_node().to_sexp());
    //while cursor.goto_first_child() {
    cursor.goto_first_child();
    // println!("{}", cursor.node().to_sexp());
    loop {
        // Advance cursor and skip if a comment
        if cursor.goto_first_child() && cursor.node().kind() != "comment" {
            // println!("{}", cursor.node().to_sexp());
            // println!("{}", cursor.node().to_sexp());

            // let cmd = if NamedCommand::try_from(cursor.node().kind()).is_ok() {
            //     // println!("{}", cursor.node().to_sexp());
            //     // TODO: add arguments
            //     CommandType::NamedCommand(NamedCommand::from_node(&cursor.node(), text)?)
            // } else {
            //     CommandType::GenericCommand(GenericCommand::from_node(&cursor.node(), text)?)
            // };
            //
            let cmd = CommandNode::from_node(&cursor.node(), text)?;

            commands.push(cmd);

            cursor.goto_parent();
        }

        // If no more commands, break!
        if !cursor.goto_next_sibling() {
            break;
        }
    }

    Ok(Ast { commands })
}

/// A command in the LAMMPS Input Script
///  TODO rename this to be a command node???
#[derive(Debug, PartialEq, Clone)]
pub struct CommandNode {
    pub command_type: CommandType,
    range: Range,
}

impl CommandNode {
    pub fn range(&self) -> Range {
        self.range
    }
}

impl FromNode for CommandNode {
    fn from_node(node: &Node, text: &[u8]) -> Result<Self, FromNodeError> {
        let range = node.range();

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
    fn from_node(node: &Node, text: &[u8]) -> Result<Self, FromNodeError> {
        let cmd = if NamedCommand::try_from(node.kind()).is_ok() {
            // println!("{}", cursor.node().to_sexp());
            // TODO: add arguments
            CommandType::NamedCommand(NamedCommand::from_node(&node, text)?)
        } else {
            CommandType::GenericCommand(GenericCommand::from_node(&node, text)?)
        };

        Ok(cmd)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct GenericCommand {
    pub name: String,
    pub args: Vec<Argument>,
    pub start: Point,
    pub end: Point,
    pub start_byte: usize,
    pub end_byte: usize,
}

impl FromNode for GenericCommand {
    fn from_node(node: &Node, text: &[u8]) -> Result<Self, FromNodeError> {
        let mut cursor = node.walk();
        // let kind = node.kind().to_string();
        let start = node.start_position();
        let end = node.end_position();
        let start_byte = node.start_byte();
        let end_byte = node.end_byte();

        let mut args = vec![];

        debug_assert!(cursor.node() == *node);

        cursor.goto_first_child();

        // TODO: use a field in the TS grammar
        let name = cursor.node().utf8_text(text)?.to_string();

        while cursor.goto_next_sibling() {
            for node in cursor.node().children(&mut cursor) {
                args.push(Argument::from_node(&node, text)?);
            }
        }

        cursor.goto_parent();
        Ok(GenericCommand {
            name,
            args,
            start,
            end,
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
    String,
    Expression(expressions::Expression),
    // TODO: Remove? Can't know if a group name until further on in the process???
    // Perhaps make it an identifier that then is decided to be either
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    UnderscoreIdent(Ident),
    Word(String),
}

impl FromNode for Argument {
    fn from_node(
        node: &Node,
        // _cursor: &mut tree_sitter::TreeCursor,
        text: &[u8],
    ) -> Result<Self, FromNodeError> {
        // TODO: make these variants more complete.
        //.child(0).unwrap()
        // Did removing above from match fix things???
        match node.kind() {
            "int" => Ok(Self::Int(node.utf8_text(text)?.parse::<isize>()?)),
            "float" => Ok(Self::Float(node.utf8_text(text)?.parse::<f64>()?)),
            "bool" => Ok(Self::Bool(match node.utf8_text(text)? {
                "on" | "yes" | "true" => true,
                "off" | "no" | "false" => false,
                name => Err(FromNodeError::UnknownCustom {
                    kind: "Bool Variant".to_owned(),
                    name: name.to_string(),
                })?,
            })),
            // TODO: Expressions not wrapped in varround are not valid???
            "expression" => Ok(Self::Expression(expressions::Expression::parse_expression(
                node, text,
            )?)),
            "string" => Ok(Self::String),
            "group" => Ok(Self::Group),
            "underscore_ident" => Ok(Self::UnderscoreIdent(Ident::new(
                &node.child(0).unwrap(),
                text,
            )?)),
            "var_curly" => Ok(Self::VarCurly(Ident::new(&node.child(1).unwrap(), text)?)),
            "var_round" => Ok(Self::VarRound(expressions::Expression::parse_expression(
                &node.child(1).unwrap(),
                text,
            )?)),
            "argname" => Ok(Self::ArgName(
                //.child(0).unwrap()
                node.utf8_text(text).unwrap().to_string(),
            )),
            "word" => Ok(Self::Word(
                //.child(0).unwrap()
                node.utf8_text(text).unwrap().to_string(),
            )),
            x => Err(FromNodeError::UnknownCustom {
                kind: "argument type".to_string(),
                name: x.to_string(),
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
            Argument::String => write!(f, "string"),
            Argument::Expression(x) => write!(f, "expression: {x}"),
            Argument::Group => write!(f, "group"),
            Argument::UnderscoreIdent(x) => write!(f, "underscore_ident: {x}"),
            Argument::Word(x) => write!(f, "word: {x}"),
        }
    }
}

/// Commands that have a special form in the tree sitter grammar
/// TODO: add arguments
/// TODO: Add command location
#[derive(Debug, PartialEq, Clone)]
pub enum NamedCommand {
    Fix(FixDef),
    Compute,
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
    fn from_node(node: &Node, text: &[u8]) -> Result<NamedCommand, FromNodeError> {
        match node.kind() {
            "fix" => Ok(NamedCommand::Fix(FixDef::from_node(node, text)?)),
            "compute" => Ok(NamedCommand::Compute),
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
            _ => Err(FromNodeError::UnknownCommand(node.kind())),
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
    pub fn range(&self) -> Range {
        self.fix_id.range()
    }
}

impl FromNode for FixDef {
    /// TODO: Hand a cursor instead???
    // TODO: Remove unwraps
    fn from_node(node: &Node, text: &[u8]) -> Result<Self, FromNodeError> {
        let mut cursor = node.walk();

        // TODO: handle case this is false
        cursor.goto_first_child();
        let mut children = node.children(&mut cursor);

        // skip the fix keyword
        children.next();
        let fix_id = Ident::new(&children.next().unwrap(), text)?;
        let group_id = children.next().unwrap().utf8_text(text)?.into();
        let fix_style = children
            .next()
            .unwrap()
            .utf8_text(text)?
            .try_into()
            .unwrap();

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

/// This doesn't work for the new case that the fix has data
impl TryFrom<&str> for NamedCommand {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "fix" => Ok(Self::Fix(FixDef::default())),
            "compute" => Ok(Self::Compute),
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
mod tests {
    use tree_sitter::Parser;

    use crate::ast::expressions::{BinaryOp, Expression};
    use crate::ast::from_node::FromNode;
    use crate::ast::{Argument, FixDef};
    use crate::fix_styles::FixStyle;
    use crate::identifinder::{Ident, IdentType};

    use super::{ts_to_ast, Ast};

    fn setup_parser() -> Parser {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        parser
    }

    #[test]
    fn test_ast() {
        let mut parser = setup_parser();
        // let source_bytes = include_bytes!("../../fix.lmp");
        let source_bytes = include_bytes!("../../in.nemd");
        let tree = parser.parse(source_bytes, None).unwrap();

        let ast = ts_to_ast(&tree, source_bytes);
        dbg!(ast.unwrap());

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
                    Argument::ArgName("temp".into()),
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
}

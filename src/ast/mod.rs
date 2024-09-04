//! Convert the treesitter trees into an AST.

// For denying unwraps and expects in this file
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]

pub mod expressions;
pub mod find_node;
pub mod from_node;

use crate::{
    spanned_error::{SpannedError, WithSpan},
    spans::{Point, Span},
};
use std::{fmt::Display, str::Utf8Error};

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

/// An AST that could not be fully parsed. Also contains errors.
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
#[derive(Debug, PartialEq, Clone)]
pub struct CommandNode {
    pub command_type: CommandType,
    range: Span,
}

impl CommandNode {
    pub fn span(&self) -> Span {
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
        // TODO: just try parsing this first.
        // TODO: Flatten the enum have generic just be the last command type.
        let cmd = if NamedCommand::try_from(node.kind()).is_ok() {
            // TODO: add arguments
            CommandType::NamedCommand(NamedCommand::from_node(node, text)?)
        } else {
            CommandType::GenericCommand(GenericCommand::from_node(node, text)?)
        };

        Ok(cmd)
    }
}

#[derive(Clone, Default, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Word {
    pub contents: String,
    pub span: Span,
}
impl Display for Word {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.contents.fmt(f)
    }
}

impl Word {
    pub fn new(contents: String, span: impl Into<Span>) -> Self {
        let contents = contents.to_owned();
        let span = span.into();

        Self { contents, span }
    }

    pub fn as_str(&self) -> &str {
        self.contents.as_str()
    }

    pub(crate) fn parse_word(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Utf8Error> {
        let contents = node.utf8_text(text.as_ref())?.to_owned();
        let span = node.range().into();

        Ok(Word { contents, span })
    }
}

impl FromNode for Word {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        Ok(Word::parse_word(node, text)?)
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct GenericCommand {
    pub name: Word,
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
        let start = node.start_position();
        let end = node.end_position();
        let start_byte = node.start_byte();
        let end_byte = node.end_byte();

        let text = text.as_ref();

        let mut args = vec![];

        debug_assert!(cursor.node() == *node);

        if !cursor.goto_first_child() {
            return Err(FromNodeError::PartialNode(
                "missing command name".to_owned(),
            ));
        }

        // TODO: use a field in the TS grammar instead?
        let name = Word::from_node(&cursor.node(), text)?;

        while cursor.goto_next_sibling() {
            for node in cursor.node().children(&mut cursor) {
                args.push(Argument::from_node(&node, text)?);
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

#[derive(Debug, PartialEq, Clone)]
pub struct Argument {
    pub kind: ArgumentKind,
    // TODO: Make the span be put on the kind rather than here?
    pub span: Span,
}

impl Argument {
    pub fn new(kind: ArgumentKind, span: impl Into<Span>) -> Argument {
        Argument {
            kind,
            span: span.into(),
        }
    }
}

/// Acceptable argument types for LAMMPS commands
#[derive(Debug, PartialEq, Clone)]
pub enum ArgumentKind {
    /// Expression from LAMMPS
    Int(isize),
    Float(f64),
    Bool(bool),
    ArgName(String), // TODO: Rename to keyword arg?
    /// Variables within curly braces
    VarCurly(Ident),
    VarRound(expressions::Expression),
    String(String),
    /// Multiple argument stuck together, usually strings and expansions.
    Concatenation(Vec<Argument>),
    /// Expression.
    /// TODO: Not really valid as a bare argument except for with variable commands
    Expression(expressions::Expression),
    // TODO: Remove? Can't know if a group name until further on in the process???
    // Perhaps make it an identifier that then is decided to be either
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    UnderscoreIdent(Ident),
    // TODO: Make the contents of this a [`Word`]
    /// A white-space delimited word
    Word(String),
    /// Indicates an invalid node in the tree-sitter grammar
    Error,
}

impl FromNode for Argument {
    type Error = FromNodeError;

    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        let text = text.as_ref();
        let kind = ArgumentKind::from_node(node, text)?;

        Ok(Argument {
            kind,
            span: node.range().into(),
        })
    }
}

impl FromNode for ArgumentKind {
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
            "concatenation" => {
                let mut cursor = node.walk();
                let v: Result<Vec<_>, _> = node
                    .children(&mut cursor)
                    .map(|ch| Argument::from_node(&ch, text))
                    .collect();
                Ok(Self::Concatenation(v?))
            }
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
        self.kind.fmt(f)
    }
}

impl Display for ArgumentKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // TODO: rework this. This is only really appropriate for debug, not display
        match self {
            Self::Int(x) => write!(f, "int: {x}"),
            Self::Float(x) => write!(f, "float: {x}"),
            Self::Bool(x) => write!(f, "bool: {x}"),
            Self::ArgName(x) => write!(f, "argname: {x}"),
            Self::VarCurly(x) => write!(f, "var_curly: ${{{0}}}", x.name),
            Self::VarRound(x) => write!(f, "var_round: $({x})"),
            Self::Concatenation(v) => {
                // TODO: this will look really ugly
                write!(f, "concatenation: ")?;
                for a in v {
                    write!(f, "{a}")?;
                }
                Ok(())
            }
            // TODO: properly implement string
            Self::String(s) => write!(f, "string: {s}"),
            Self::Expression(x) => write!(f, "expression: {x}"),
            Self::Group => write!(f, "group"),
            Self::UnderscoreIdent(x) => write!(f, "underscore_ident: {x}"),
            Self::Word(x) => write!(f, "word: {x}"),
            Self::Error => write!(f, "ERROR"),
        }
    }
}

/// Commands that have a special form in the tree sitter grammar
/// TODO: add arguments
/// TODO: Add command location
/// TODO: Reduce the number of these, both here and in the grammar.
#[derive(Debug, PartialEq, Clone)]
pub enum NamedCommand {
    /// A Fix defintion
    Fix(FixDef),
    /// A compute definition
    Compute(ComputeDef),
    Style,
    Modify,
    AtomStyle,
    Boundary,
    /// A variable definition
    VariableDef(VariableDef),
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
            "variable_def" => Ok(NamedCommand::VariableDef(VariableDef::from_node(
                node, text,
            )?)),
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
    pub fix_id: Ident,  // Or just keep as a string?
    pub group_id: Word, // TODO:  Create group identifiers
    pub fix_style: FixStyle,
    // Arguments for fix command
    pub args: Vec<Argument>,
    pub span: Span,
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
        let span = node.range().into();
        let mut cursor = node.walk();
        let text = text.as_ref();

        let mut children = node.children(&mut cursor);

        // skip the fix keyword
        children.next();
        let fix_id = Ident::new(
            &children
                .next()
                .ok_or(FromNodeError::PartialNode("Missing fix identifier".into()))?,
            text,
        )?;

        let group_id = Word::from_node(&children.next().into_err()?, text)?;
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
            span,
        })
    }
}

#[derive(Debug, Default, PartialEq, Clone)]
pub struct ComputeDef {
    pub compute_id: Ident, // Or just keep as a string?
    pub group_id: Word,    // TODO:  Create group identifiers
    pub compute_style: ComputeStyle,
    // Arguments for fix command
    pub args: Vec<Argument>,
    pub span: Span,
}

impl ComputeDef {
    pub fn range(&self) -> Span {
        self.span
    }
}

impl FromNode for ComputeDef {
    type Error = FromNodeError;
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        // Note: Might be more efficient to hand a cursor instead.
        let span = node.range().into();
        let mut cursor = node.walk();

        let mut children = node.children(&mut cursor);
        let text = text.as_ref();

        // skip the fix keyword
        children.next();
        let compute_id = Ident::new(&children.next().into_err()?, text)?;
        let group_id = Word::from_node(&children.next().into_err()?, text)?;
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
            span,
        })
    }
}

//TODO: finish me
#[derive(Debug, Default, PartialEq, Clone)]
pub struct VariableDef {
    pub variable_id: Ident, // Or just keep as a string?
    pub variable_style: Word,
    // Arguments for fix command
    pub args: Vec<Argument>,
    pub span: Span,
}

impl FromNode for VariableDef {
    type Error = FromNodeError;

    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        let text = text.as_ref();
        let span: Span = node.range().into();
        let mut cursor = node.walk();

        let mut children = node.children(&mut cursor);

        let variable_keyword = children
            .next()
            // Should basically be impossible
            .ok_or(Self::Error::PartialNode(
                "missing variable keyword".to_string(),
            ))?;

        debug_assert_eq!(variable_keyword.utf8_text(text)?, "variable");

        let variable_ident = children
            .next()
            // Should basically be impossible
            .ok_or(Self::Error::PartialNode(
                "missing variable identifier".to_string(),
            ))?;

        let variable_ident = dbg!(Ident::new(&variable_ident, text)?);

        // TODO: Make this missing thing nicer.
        let variable_kind = Word::from_node(
            &children.next().ok_or(Self::Error::PartialNode(
                "missing variable style".to_string(),
            ))?,
            text,
        )?;

        let args: Result<Vec<Argument>, _> = children
            .map(|arg| Argument::from_node(&arg, text))
            .collect();

        let args = args?;

        if args.is_empty() {
            return Err(Self::Error::PartialNode(
                "missing arguments in variable command".to_string(),
            ));
        }

        // TODO: Ensure that the style is valid.

        // These styles expect just a single expression.
        if matches!(variable_kind.contents.as_str(), "equal" | "vector" | "atom") {
            match args[0].kind {
                ArgumentKind::Expression(_) => (),
                _ => {
                    return Err(Self::Error::PartialNode(format!(
                        "expected expression for variable style {}",
                        variable_kind.contents
                    )));
                }
            }

            if args.len() > 1 {
                return Err(Self::Error::PartialNode(format!(
                    "only one argument expected for variable style {}",
                    variable_kind.contents
                )));
            }
        }

        Ok(VariableDef {
            variable_id: variable_ident,
            variable_style: variable_kind,
            args,
            span,
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
            "variable_def" => Ok(Self::VariableDef(VariableDef::default())),
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

    use core::panic;

    use pretty_assertions::assert_eq;
    use tree_sitter::Parser;

    use crate::ast::expressions::{BinaryOp, Expression};
    use crate::ast::from_node::FromNode;
    use crate::ast::{Argument, ComputeDef, FixDef, NamedCommand, VariableDef};
    use crate::ast::{ArgumentKind, Word};
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

        let _ast = ts_to_ast(&tree, source_bytes);
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
                    span: ((0, 4)..(0, 7)).into()
                },
                group_id: Word {
                    contents: "all".into(),
                    span: ((0, 8)..(0, 11)).into()
                },
                fix_style: FixStyle::Nve,
                args: vec![],
                span: ((0, 0)..(0, 15)).into(),
            }
        );
    }

    #[test]
    // TODO: Finish this test
    fn parse_index_variable() {
        let mut parser = setup_parser();
        let text = b"variable file_name index step4.1.atm\n";

        let tree = parser.parse(text, None).unwrap();

        dbg!(tree.root_node().to_sexp());
        let cmd_node = tree.root_node().child(0).unwrap(); // command node
        let var_def_node = cmd_node.child(0).unwrap(); // variable node

        dbg!(tree.root_node().to_sexp());

        let expected = VariableDef {
            variable_id: Ident {
                name: "file_name".to_string(),
                ident_type: IdentType::Variable,
                span: ((0, 9)..(0, 18)).into(),
            },
            variable_style: Word::new("index".to_string(), (0, 19)..(0, 24)),
            args: vec![Argument {
                kind: ArgumentKind::Word("step4.1.atm".to_string()),
                span: ((0, 25)..(0, 36)).into(),
            }],
            span: ((0, 0)..(0, 36)).into(),
        };

        assert_eq!(
            VariableDef::from_node(&var_def_node, text).as_ref(),
            Ok(&expected)
        );

        let ast = ts_to_ast(&tree, text);
        assert!(ast.is_ok());

        if let Ok(ast) = ast {
            assert_eq!(ast.commands.len(), 1);
            match &ast.commands[0].command_type {
                crate::ast::CommandType::NamedCommand(NamedCommand::VariableDef(var)) => {
                    assert_eq!(*var, expected)
                }
                cmd => panic!("Unexpected command {cmd:?}"),
            }
        }
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
        use ArgumentKind as AK;
        assert_eq!(
            // ast.commands[0],
            FixDef::from_node(&fix_node, source_bytes.as_slice()).unwrap(),
            FixDef {
                fix_id: Ident {
                    name: "NVT".into(),
                    ident_type: IdentType::Fix,
                    span: ((0, 4)..(0, 7)).into()
                },
                group_id: Word {
                    contents: "all".into(),
                    span: ((0, 8)..(0, 11)).into()
                },
                fix_style: FixStyle::Nvt,
                args: vec![
                    Argument::new(AK::Word("temp".into()), (0, 16)..(0, 20)),
                    Argument::new(AK::Int(1), (0, 21)..(0, 22)),
                    Argument::new(AK::Float(1.5), (0, 23)..(0, 26)),
                    Argument::new(
                        AK::VarRound(Expression::BinaryOp(
                            Expression::Float(100.0).into(),
                            BinaryOp::Multiply,
                            Expression::ThermoKeyword(Word {
                                contents: "dt".into(),
                                span: ((0, 35)..(0, 37)).into()
                            })
                            .into(),
                        )),
                        (0, 27)..(0, 38)
                    )
                ],
                span: ((0, 0)..(0, 38)).into()
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
                    span: ((0, 8), (0, 13)).into()
                },
                group_id: Word {
                    contents: "water".into(),
                    span: ((0, 14), (0, 19)).into()
                },
                compute_style: ComputeStyle::TempRegion,
                // TODO: Change to a more generic word argument
                args: vec![Argument::new(
                    ArgumentKind::Word("hot_region".into()),
                    ((0, 32), (0, 42))
                ),],
                span: ((0, 0), (0, 42)).into()
            }
        );
    }

    #[test]
    fn parse_variable_def() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 1.0*dt";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();

        let expected = VariableDef {
            variable_id: Ident {
                name: "a".into(),
                ident_type: IdentType::Variable,
                // TODO: Check this is the type expected.
                span: ((0, 9), (0, 10)).into(),
            },
            variable_style: Word {
                contents: "equal".into(),
                span: ((0, 11), (0, 16)).into(),
            },

            args: vec![Argument {
                kind: ArgumentKind::Expression(Expression::BinaryOp(
                    Box::new(Expression::Float(1.0)),
                    BinaryOp::Multiply,
                    Box::new(Expression::ThermoKeyword(Word {
                        contents: "dt".to_owned(),
                        span: ((0, 21), (0, 23)).into(),
                    })),
                )),
                span: ((0, 17), (0, 23)).into(),
            }],
            span: ((0, 0), (0, 23)).into(),
        };

        let variable_node = command_node.child(0).expect("Should find child node.");
        let parsed = VariableDef::from_node(&variable_node, source_bytes);
        assert_eq!(parsed, Ok(expected));
    }
}

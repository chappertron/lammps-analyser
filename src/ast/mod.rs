//! Convert the treesitter trees into an AST.

// For denying unwraps and expects in this file
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]

pub mod expressions;
pub mod find_node;
pub mod from_node;
pub mod impl_helpers;

use crate::{
    spanned_error::SpannedError,
    spans::{Point, Span},
};
use std::{fmt::Display, str::Utf8Error};

use expressions::Expression;
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

/// TODO: return both the partial AST and the Errors
pub fn ts_to_ast(tree: &Tree, text: impl AsRef<[u8]>) -> Result<Ast, PartialAst> {
    let mut cursor = tree.walk();
    let text = text.as_ref();

    let mut commands = vec![];
    let mut errors = vec![];
    cursor.goto_first_child();

    for node in tree.root_node().named_children(&mut cursor) {
        // TODO: This if-statement may be removeable
        if node.kind() != "comment" {
            let result = CommandNode::from_node(&node, text); //.with_span(node.range().into());

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
    type Error = (Self, FromNodeError);
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, Self::Error> {
        let range = node.range().into();

        match CommandType::from_node(node, text) {
            Ok(command_type) => Ok(CommandNode {
                command_type,
                range,
            }),
            Err((command_type, err)) => Err((
                CommandNode {
                    command_type,
                    range,
                },
                err,
            )),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub enum CommandType {
    GenericCommand(GenericCommand),
    /// A Fix defintion
    Fix(FixDef),
    /// A compute definition
    Compute(ComputeDef),
    /// A variable definition
    VariableDef(VariableDef),
    Shell,
    Error,
}

impl FromNode for CommandType {
    // NOTE: This is chosen as the error type so an error node can be returned instead
    type Error = (Self, FromNodeError);
    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> Result<Self, (Self, FromNodeError)> {
        let mut result = match node.kind() {
            "fix" => Ok(Self::Fix(
                FixDef::from_node(node, &text).map_err(|e| (Self::Error, e))?,
            )),
            "compute" => Ok(Self::Compute(
                ComputeDef::from_node(node, &text).map_err(|e| (Self::Error, e))?,
            )),
            "variable_def" => Ok(Self::VariableDef(
                VariableDef::from_node(node, &text).map_err(|e| (Self::Error, e))?,
            )),
            "shell" => Ok(Self::Shell),
            // Fall back to the generic command type
            "command" => Ok(Self::GenericCommand(
                GenericCommand::from_node(node, &text).map_err(|e| (Self::Error, e))?,
            )),
            _ => Ok(Self::Error),
        };

        if let Ok(Self::Error) = result {
            // NOTE: Check the child node and infer what type of node it was supposed to be
            // and pass to that parser.

            result = match node.child(0).map(|node| node.kind()) {
                // Try and pass this as a compute
                Some("compute") => Ok(Self::Compute(
                    ComputeDef::from_node(node, &text).map_err(|e| (Self::Error, e))?,
                )),
                Some("fix_id") => Ok(Self::Fix(
                    FixDef::from_node(node, &text).map_err(|e| (Self::Error, e))?,
                )),
                Some("variable") => Ok(Self::VariableDef(
                    VariableDef::from_node(node, &text).map_err(|e| (Self::Error, e))?,
                )),
                // Cannot further process
                // NOTE: This is `Ok` rather than `Error` because syntax errors are also
                // found within `ErrorFinder`
                _ => Ok(Self::Error),
            };
        }

        result
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
    /// TODO: Make a quoted expression a seperate thing?
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
            "quoted_expression" => quoted_expression(node, text).map(Self::Expression),
            "string" => Ok(Self::String(node.utf8_text(text)?.to_owned())),
            "group" => Ok(Self::Group),
            "underscore_ident" => Ok(Self::UnderscoreIdent(Ident::new(
                &node.child(0).into_err()?,
                text,
            )?)),
            // TODO: check surrounding brackets
            "var_curly" => Ok(Self::VarCurly(Ident::new(
                &node.child(1).into_err()?,
                text,
            )?)),
            // TODO: check surrounding brackets
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

// TODO: Move this to a more sensible module
pub(crate) fn quoted_expression(node: &Node, text: &[u8]) -> Result<Expression, FromNodeError> {
    let mut cursor = node.walk();
    let mut children = node.children(&mut cursor);

    let is_left_missing = match children.next() {
        Some(left_quote) if left_quote.kind() != "\"" => true,
        None => true,
        Some(_) => false,
    };

    if is_left_missing {
        return Err(FromNodeError::PartialNode("expected '\"'".into()));
    }

    let expr = children
        .next()
        .ok_or(FromNodeError::PartialNode("expected expression".into()))?;

    let expr = Expression::parse_expression(&expr, text)?;

    let is_right_missing = match children.next() {
        Some(right_quote) if right_quote.kind() != "\"" => true,
        None => true,
        Some(_) => false,
    };

    if is_right_missing {
        return Err(FromNodeError::PartialNode("expected '\"'".into()));
    }

    Ok(expr)
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

        // skip the compute keyword
        children.next();
        let compute_id = children
            .next()
            .ok_or(Self::Error::PartialNode("missing compute ID".to_string()))?;

        let compute_id = Ident::new(&compute_id, text)?;

        let group_id = Word::from_node(
            &children
                .next()
                .ok_or(Self::Error::PartialNode("missing group ID".to_string()))?,
            text,
        )?;
        let compute_style = children
            .next()
            .ok_or(Self::Error::PartialNode(
                "missing compute style".to_string(),
            ))?
            .utf8_text(text)?
            .into();

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

        let variable_ident = Ident::new(&variable_ident, text)?;

        // TODO: Make this missing thing nicer.
        let variable_kind = Word::from_node(
            &children.next().ok_or(Self::Error::PartialNode(
                "missing variable style".to_string(),
            ))?,
            text,
        )?;

        // TODO: Check if is missing node here?
        let args: Result<Vec<Argument>, _> = children
            .map(|arg| {
                if arg.is_missing() {
                    return Err(FromNodeError::PartialNode(
                        "missing variable expression".to_string(),
                    ));
                }
                Argument::from_node(&arg, text)
            })
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

#[cfg(test)]
// Allow unwraps in the tests module, but not in the parent module.
#[allow(clippy::unwrap_used)]
#[allow(clippy::expect_used)]
mod tests {

    use pretty_assertions::assert_eq;

    use crate::ast::expressions::{BinaryOp, Expression};
    use crate::ast::from_node::{FromNode, FromNodeError};
    use crate::ast::{Argument, ComputeDef, FixDef, VariableDef};
    use crate::ast::{ArgumentKind, Word};
    use crate::compute_styles::ComputeStyle;
    use crate::fix_styles::FixStyle;
    use crate::identifinder::{Ident, IdentType};
    use crate::utils::testing::parse;

    use super::ts_to_ast;

    #[test]
    #[ignore = "Not yet fully implemented"]
    fn test_ast() {
        // let source_bytes = include_bytes!("../../fix.lmp");
        let source_bytes = include_bytes!("../../example_input_scripts/in.nemd");
        let tree = parse(source_bytes);

        let _ast = ts_to_ast(&tree, source_bytes);
        // dbg!(ast.unwrap());

        unimplemented!()
    }

    #[test]
    fn parse_fix_no_args() {
        let source_bytes = b"fix NVE all nve";

        let tree = parse(source_bytes);

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        // Lots of tedium to parsing this...
        let fix_node = dbg!(root_node.child(0)).unwrap();

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
    fn parse_index_variable() {
        let text = "variable file_name index step4.1.atm";
        // let text = include_str!("../../example_input_scripts/in.variable_index");

        let tree = parse(text);

        dbg!(tree.root_node().to_sexp());
        let root_node = tree.root_node();
        let var_def_node = root_node.child(0).unwrap(); // variable node

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
            dbg!(&ast);
            assert_eq!(ast.commands.len(), 1);
            match &ast.commands[0].command_type {
                crate::ast::CommandType::VariableDef(var) => {
                    assert_eq!(*var, expected)
                }
                cmd => panic!("Unexpected command {cmd:?}"),
            }
        }
    }

    #[test]
    fn parse_fix_with_args() {
        let source_bytes = b"fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = parse(source_bytes);

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let fix_node = dbg!(root_node.child(0)).unwrap();

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
        let source_bytes = b"compute T_hot water temp/region hot_region";

        let tree = parse(source_bytes);

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let compute_node = dbg!(root_node.child(0)).unwrap();

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
        let source_bytes = b"variable a equal 1.0*dt";

        let tree = parse(source_bytes);

        let root_node = tree.root_node();

        let expected = VariableDef {
            variable_id: Ident {
                name: "a".into(),
                ident_type: IdentType::Variable,
                // TODO: Check this is the type expected.
                span: ((0, 9)..(0, 10)).into(),
            },
            variable_style: Word {
                contents: "equal".into(),
                span: ((0, 11)..(0, 16)).into(),
            },

            args: vec![Argument {
                kind: ArgumentKind::Expression(Expression::BinaryOp(
                    Box::new(Expression::Float(1.0)),
                    BinaryOp::Multiply,
                    Box::new(Expression::ThermoKeyword(Word {
                        contents: "dt".to_owned(),
                        span: ((0, 21)..(0, 23)).into(),
                    })),
                )),
                span: ((0, 17)..(0, 23)).into(),
            }],
            span: ((0, 0)..(0, 23)).into(),
        };

        let variable_node = root_node.child(0).expect("Should find child node.");
        let parsed = VariableDef::from_node(&variable_node, source_bytes);
        assert_eq!(parsed, Ok(expected));
    }

    #[test]
    fn parse_incomplete_variable_def() {
        let source_bytes = b"variable a equal";

        let tree = parse(source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        assert!(!root_node.is_missing());

        let expected =
            FromNodeError::PartialNode("missing arguments in variable command".to_string());

        let variable_node = root_node.child(0).expect("Should find child node.");
        let parsed = VariableDef::from_node(&variable_node, source_bytes);
        assert_eq!(parsed, Err(expected));
    }

    #[test]
    fn parse_incomplete_compute_def() {
        let source_bytes = b"compute a ";

        let tree = parse(source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        assert!(!root_node.is_missing());

        let ast = ts_to_ast(&tree, source_bytes);

        // TODO: double check more about the syntax tree.
        assert!(ast.is_err());

        let compute_node = root_node.child(0).expect("Should find child node.");
        let parsed = ComputeDef::from_node(&compute_node, source_bytes);
        use crate::ast::from_node::FromNodeError;
        assert_eq!(
            parsed,
            Err(FromNodeError::PartialNode("missing group ID".to_string()))
        );
    }
}

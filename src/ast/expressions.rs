// Deny unwraps and expects in this file to prevent crashes in the LSP
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]

use std::fmt::Display;

use tree_sitter::Node;

use crate::{identifinder::Ident, utils::into_error::IntoError};
use std::convert::TryFrom;
use thiserror::Error;

use super::{from_node::MissingNode, Word};

#[derive(Debug, Default, PartialEq, Clone)]
/// An mathematical expression in LAMMPS
pub enum Expression {
    #[default]
    Missing,
    /// LAMMPS Identifier for a fix/compute/variable that is
    /// proceeded by f_/c_/v_ to indicate the type.
    UnderscoreIdent(Ident),
    /// An integer
    Int(isize),
    /// A float. Parsed as a 64-bit float
    Float(f64),
    /// Boolean Operation
    Bool(bool),
    /// Unary operation
    UnaryOp(UnaryOp, Box<Expression>),
    /// A binary expression between two other expressions.
    BinaryOp(Box<Expression>, BinaryOp, Box<Expression>),
    /// A built-in function. Slightly deviates from grammar, by only having one function type.
    // TODO: Make a cow for name?
    Function(Word, Vec<Expression>),
    /// An expression wrapped in brackets.
    Parens(Box<Expression>),
    /// Thermo keyword
    /// TODO: Use an enum instead
    ThermoKeyword(Word),
    /// TODO: Word might not actually be implemented for expr
    Word(Word),
}

impl Display for Expression {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Missing => write!(f, ""),
            Self::UnderscoreIdent(ident) => write!(f, "{}", ident.underscore_ident()),
            Self::Int(i) => write!(f, "{i}"),
            Self::Float(fl) => write!(f, "{fl}"),
            Self::Bool(b) => write!(f, "{b}"),
            Self::UnaryOp(op, expr) => write!(f, "{op}{expr}"),
            Self::BinaryOp(lhs, op, rhs) => write!(f, "{lhs} {op} {rhs}"),
            Self::Function(name, args) => {
                write!(f, "{name}")?;
                write!(f, "(")?;
                for (i, arg) in args.iter().enumerate() {
                    write!(f, "{}", arg)?;
                    if i != args.len() - 1 {
                        write!(f, ", ")?;
                    }
                }
                write!(f, ")")
            }
            Self::Parens(expr) => write!(f, "({expr})"),
            Self::ThermoKeyword(kw) => write!(f, "{kw}"),
            Self::Word(word) => write!(f, "{word}"),
        }
    }
}

// TODO: Flatten this into the FromNodeError?
// This means that conversion from word can't be done
#[allow(clippy::enum_variant_names)]
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ParseExprError {
    #[error("Could not parse text as UTF-8 {0}")]
    Utf8Error(std::str::Utf8Error),
    #[error("Could not parse text as int {0}")]
    ParseIntError(std::num::ParseIntError),
    #[error("Could not parse text as float {0}")]
    ParseFloatError(std::num::ParseFloatError),
    #[error("Invalid unary operator: {0}")]
    InvalidUnaryOperator(String), // TODO: make &'a str
    #[error("Invalid binary operator: {0}")]
    InvalidBinaryOperator(String), // TODO: Make &'a str
    #[error("Tree-sitter Error Node Found")]
    ErrorNode,
    #[error("Unknown Expression type: {0}")]
    UnknownExpressionType(String),
    #[error("Syntax Error: expected expression")]
    MissingToken,
}

impl From<std::num::ParseFloatError> for ParseExprError {
    fn from(v: std::num::ParseFloatError) -> Self {
        Self::ParseFloatError(v)
    }
}

impl From<std::num::ParseIntError> for ParseExprError {
    fn from(v: std::num::ParseIntError) -> Self {
        Self::ParseIntError(v)
    }
}

impl From<std::str::Utf8Error> for ParseExprError {
    fn from(v: std::str::Utf8Error) -> Self {
        Self::Utf8Error(v)
    }
}
impl From<MissingNode> for ParseExprError {
    fn from(_: MissingNode) -> Self {
        // TODO: Add more context here. What went wrong?
        Self::MissingToken
    }
}

impl Expression {
    /// TODO: Handle Errors
    pub(crate) fn parse_expression(node: &Node<'_>, text: &[u8]) -> Result<Self, ParseExprError> {
        // TODO: Handle missing node
        if node.is_missing() {
            return Err(ParseExprError::MissingToken);
        }

        if let Some(child) = node.child(0) {
            if child.is_missing() {
                return Err(ParseExprError::MissingToken);
            }
        }

        match node.kind() {
            // Child 0 is the identifier
            // The v_/f_/c_ prefix is anonymous
            "underscore_ident" => Ok(Self::UnderscoreIdent(Ident::new(
                &node.child(0).into_err()?,
                text,
            )?)),
            "binary_op" => Ok(Self::BinaryOp(
                Box::new(Self::parse_expression(&node.child(0).into_err()?, text)?),
                // TODO: Find a way to get the operator from the TS symbol
                BinaryOp::try_from(node.child(1).into_err()?.utf8_text(text)?)?,
                Box::new(Self::parse_expression(&node.child(2).into_err()?, text)?),
            )),

            "unary_op" => Ok(Self::UnaryOp(
                // TODO: Can the operator be parsed with the kind??
                UnaryOp::try_from(node.child(0).into_err()?.utf8_text(text)?)?,
                Self::parse_expression(&node.child(1).into_err()?, text)?.into(),
            )),
            "binary_func" | "unary_func" | "ternary_func" | "hexnary_func" | "other_func"
            | "group_func" | "region_func" => {
                let mut cursor = node.walk();

                let mut args = Vec::new();
                let name_node = &node.child(0).into_err()?;
                let name = Word::parse_word(name_node, text)?;

                // child 0 = function name
                // child 1 = opening bracket
                for node in node.children(&mut cursor).skip(2) {
                    node.utf8_text(text)?;
                    if node.kind() == ")" {
                        break;
                    }
                    if node.kind() == "," {
                        continue;
                    }
                    args.push(Self::parse_expression(&node, text)?);
                }
                Ok(Self::Function(name, args))
            }

            "int" => Ok(Self::Int(node.utf8_text(text)?.parse()?)),

            "float" => Ok(Self::Float(node.utf8_text(text)?.parse()?)),

            // Just go into next level down
            "expression" => Ok(Self::parse_expression(&node.child(0).into_err()?, text)?),
            "parens" => Ok(Self::Parens(Box::new(Self::parse_expression(
                &node.child(1).into_err()?,
                text,
            )?))),
            "bool" => match node.utf8_text(text)? {
                "true" | "on" | "yes" => Ok(Self::Bool(true)),
                "false" | "off" | "no" => Ok(Self::Bool(false)),
                _ => unreachable!(), // TODO: Verify this is the case?
            },
            "thermo_kwarg" => Ok(Self::ThermoKeyword(Word::parse_word(node, text)?)),
            "word" => Ok(Self::Word(Word::parse_word(node, text)?)),
            "ERROR" => Err(ParseExprError::ErrorNode),
            x => Err(ParseExprError::UnknownExpressionType(x.to_owned())),
        }
        // todo!()
    }
}

#[derive(Debug, Default, PartialEq, Eq, Clone, Copy)]
pub enum BinaryOp {
    #[default]
    Add,
    Subtract,
    Multiply,
    Divide,
    Power,
    Modulo,
    Equal,
    NotEqual,
    LessThan,
    LessThanOrEqual,
    GreaterThan,
    GreaterThanOrEqual,
    And,
    Or,
    Xor,
}

impl TryFrom<&str> for BinaryOp {
    type Error = ParseExprError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "+" => Ok(Self::Add),
            "-" => Ok(Self::Subtract),
            "*" => Ok(Self::Multiply),
            "/" => Ok(Self::Divide),
            "^" => Ok(Self::Power),
            "%" => Ok(Self::Modulo),
            "==" => Ok(Self::Equal),
            "!=" => Ok(Self::NotEqual),
            "<" => Ok(Self::LessThan),
            "<=" => Ok(Self::LessThanOrEqual),
            ">" => Ok(Self::GreaterThan),
            ">=" => Ok(Self::GreaterThanOrEqual),
            "&&" => Ok(Self::And),
            "||" => Ok(Self::Or),
            "^|" => Ok(Self::Xor),
            _ => Err(ParseExprError::InvalidBinaryOperator(value.to_owned())),
        }
    }
}

impl Display for BinaryOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Add => write!(f, "+"),
            Self::Subtract => write!(f, "-"),
            Self::Multiply => write!(f, "*"),
            Self::Divide => write!(f, "/"),
            Self::Power => write!(f, "^"),
            Self::Modulo => write!(f, "%"),
            Self::Equal => write!(f, "=="),
            Self::NotEqual => write!(f, "!="),
            Self::LessThan => write!(f, "<"),
            Self::LessThanOrEqual => write!(f, "<="),
            Self::GreaterThan => write!(f, ">"),
            Self::GreaterThanOrEqual => write!(f, ">="),
            Self::And => write!(f, "&&"),
            Self::Or => write!(f, "||"),
            Self::Xor => write!(f, "^|"),
        }
    }
}

#[derive(Debug, Default, PartialEq, Eq, Clone, Copy)]
pub enum UnaryOp {
    #[default]
    Negate,
    Not,
}
impl Display for UnaryOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            UnaryOp::Negate => write!(f, "-"),
            UnaryOp::Not => write!(f, "!"),
        }
    }
}

impl TryFrom<&str> for UnaryOp {
    type Error = ParseExprError;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "-" => Ok(Self::Negate),
            "!" => Ok(Self::Not),
            value => Err(Self::Error::InvalidUnaryOperator(value.to_owned())),
        }
    }
}

#[cfg(test)]
// Allow unwraps in the tests module, but not in the parent module.
#[allow(clippy::unwrap_used)]
#[allow(clippy::expect_used)]
mod tests {

    use pretty_assertions::assert_eq;

    use crate::identifinder::IdentType;

    use super::*;

    use tree_sitter::Parser;
    fn setup_parser() -> Parser {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        parser
    }
    #[test]
    fn parse_expr() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 1+2";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Int(1)),
                BinaryOp::Add,
                Box::new(Expression::Int(2))
            )
        );
    }

    #[test]
    fn parse_expr_floats() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 1.0+2.0";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Float(1.0)),
                BinaryOp::Add,
                Box::new(Expression::Float(2.0))
            )
        );
    }

    #[test]
    fn parse_nested_expr() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 0.5*(1.0+2.0)";

        let tree = parser.parse(source_bytes, None).unwrap();

        dbg!(tree.root_node().to_sexp());
        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Float(0.5)),
                BinaryOp::Multiply,
                Box::new(Expression::Parens(Box::new(Expression::BinaryOp(
                    Box::new(Expression::Float(1.0)),
                    BinaryOp::Add,
                    Box::new(Expression::Float(2.0))
                ))))
            )
        );
    }

    #[test]
    fn parse_expr_vars() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 0.5*(v_example+temp)";

        let tree = parser.parse(source_bytes, None).unwrap();

        dbg!(tree.root_node().to_sexp());
        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Float(0.5)),
                BinaryOp::Multiply,
                Box::new(Expression::Parens(Box::new(Expression::BinaryOp(
                    Box::new(Expression::UnderscoreIdent(Ident {
                        name: "example".into(),
                        ident_type: IdentType::Variable,
                        span: ((0, 24)..(0, 31)).into()
                    })),
                    BinaryOp::Add,
                    Box::new(Expression::ThermoKeyword(Word {
                        contents: "temp".to_owned(),
                        span: ((0, 32), (0, 36)).into()
                    }))
                ))))
            )
        );
    }

    #[test]
    fn parse_unary_func() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal abs(v_example)";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()),
            Ok(Expression::Function(
                Word {
                    contents: "abs".into(),
                    span: ((0, 17), (0, 20)).into()
                },
                vec!(Expression::UnderscoreIdent(Ident {
                    name: "example".into(),
                    ident_type: IdentType::Variable,
                    span: ((0, 23), (0, 30)).into()
                })),
            ))
        );
    }

    #[test]
    fn parse_binary_func() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal ramp(v_example, 3.0)";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::Function(
                Word {
                    contents: "ramp".into(),
                    span: ((0, 17), (0, 21)).into()
                },
                vec![
                    Expression::UnderscoreIdent(Ident {
                        name: "example".into(),
                        ident_type: IdentType::Variable,
                        span: ((0, 22), (0, 31)).into()
                    }),
                    Expression::Float(3.0)
                ],
            )
        );
    }

    #[test]
    fn func_display() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal ramp(v_example, 3.0)";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(root_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        let func = Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap();
        assert_eq!(
            // ast.commands[0],
            func,
            Expression::Function(
                Word {
                    contents: "ramp".into(),
                    span: ((0, 17)..(0, 21)).into()
                },
                vec![
                    Expression::UnderscoreIdent(Ident {
                        name: "example".into(),
                        ident_type: IdentType::Variable,
                        span: ((0, 24)..(0, 31)).into()
                    }),
                    Expression::Float(3.0)
                ],
            )
        );
        assert_eq!(func.to_string(), "ramp(v_example, 3)")
    }
}

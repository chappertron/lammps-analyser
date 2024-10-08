use std::fmt::Display;

use tree_sitter::Node;

use crate::{
    ast::Ident,
    spans::Span,
    utils::{
        into_error::IntoError,
        tree_sitter_helpers::{ExpectNode, NodeExt},
    },
};

use super::{
    expressions::Index,
    from_node::{FromNode, FromNodeError},
    Expression,
};

#[derive(Debug, PartialEq, Clone)]
pub struct Argument {
    pub kind: ArgumentKind,
    // TODO: Make the span be put on the kind rather than here?
    pub span: Span,
}

/// Acceptable argument types for LAMMPS commands
#[derive(Debug, PartialEq, Clone)]
pub enum ArgumentKind {
    /// Expression from LAMMPS
    Int(isize),
    Float(f64),
    Bool(bool),
    ArgName(String), // TODO: Is this still in use in the tree-sitter grammar?
    /// Variable expansion within curly braces
    VarCurly(Ident),
    /// A simple variable expansion in $var
    SimpleExpansion(Ident),
    /// Expression evaluations in $()
    VarRound(Expression),
    String(String),
    /// Multiple argument stuck together, usually strings and expansions.
    Concatenation(Vec<Argument>),
    /// Expression.
    /// TODO: Not really valid as a bare argument except for with variable commands
    /// TODO: Make a quoted expression a separate thing for validation.
    Expression(Expression),
    // TODO: Remove? Can't know if a group name until further on in the process???
    // Perhaps make it an identifier that then is decided to be either
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    UnderscoreIdent(Ident),
    /// An indexed `UnderscoreIdent`
    IndexedIdent(Ident, Index),

    // TODO: Make the contents of this a [`Word`]
    /// A white-space delimited word
    Word(String),
    /// Indicates an invalid node in the tree-sitter grammar
    Error,
}

impl Argument {
    pub fn new(kind: ArgumentKind, span: impl Into<Span>) -> Argument {
        Argument {
            kind,
            span: span.into(),
        }
    }
}

impl FromNode for Argument {
    type Error = FromNodeError;

    fn from_node(node: &Node, text: &str) -> Result<Self, Self::Error> {
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
        text: &str,
    ) -> Result<Self, <Argument as FromNode>::Error> {
        // TODO: make these variants more complete.
        let text = text.as_ref();
        match node.kind() {
            "int" => Ok(Self::Int(node.str_text(text).parse::<isize>()?)),
            "float" => Ok(Self::Float(node.str_text(text).parse::<f64>()?)),
            "bool" => Ok(Self::Bool(match node.str_text(text) {
                "on" | "yes" | "true" => true,
                "off" | "no" | "false" => false,
                invalid_name => Err(FromNodeError::UnknownCustom {
                    kind: "Bool Variant".to_owned(),
                    name: invalid_name.to_string(),
                    start: node.start_position(),
                })?,
            })),
            // TODO: Expressions not wrapped in varround are not valid? Add a lint.
            "expression" => Ok(Self::Expression(Expression::parse_expression(node, text)?)),
            "quoted_expression" => quoted_expression(node, text).map(Self::Expression),
            "string" => Ok(Self::String(node.str_text(text).to_owned())),
            "group" => Ok(Self::Group),
            "underscore_ident" => Ok(Self::UnderscoreIdent(underscore_ident(node, text)?)),
            "simple_expansion" => Ok(Self::SimpleExpansion(Ident::new(
                &node.child(1).into_err()?,
                text,
            )?)),
            // TODO: check surrounding brackets
            "var_curly" => Ok(Self::VarCurly(Ident::new(
                &node.child(1).into_err()?,
                text,
            )?)),
            // TODO: check surrounding brackets
            "var_round" => Ok(Self::VarRound(Expression::parse_expression(
                &node.child(1).into_err()?,
                text,
            )?)),
            "argname" => Ok(Self::ArgName(
                //.child(0).into_err()?
                node.str_text(text).to_string(),
            )),
            "word" => Ok(Self::Word(
                //.child(0).into_err()?
                node.str_text(text).to_string(),
            )),
            "concatenation" => {
                let mut cursor = node.walk();
                let v: Result<Vec<_>, _> = node
                    .children(&mut cursor)
                    .map(|ch| Argument::from_node(&ch, text))
                    .collect();
                Ok(Self::Concatenation(v?))
            }
            "indexed_ident" => {
                let mut cursor = node.walk();
                let mut children = node.children(&mut cursor);

                let ident = children
                    .next()
                    .expect_node("expected underscore identifier")?;
                let ident = underscore_ident(&ident, text)?;

                _ = children.next().expect_kind("[", "expected `[`")?;

                let index = children
                    .next()
                    .expect_node("expected `*` or integer index")?;

                let index = Index::from_node(&index, text)?;

                _ = children.next().expect_kind("]", "expected `]`")?;

                Ok(ArgumentKind::IndexedIdent(ident, index))
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

fn underscore_ident(node: &Node, text: &str) -> Result<Ident, FromNodeError> {
    // 0 is v_
    // 1 is x
    Ok(Ident::new(&node.child(1).into_err()?, text)?)
}

/// Extract a quoted expressison from a node
fn quoted_expression(node: &Node, text: &str) -> Result<Expression, FromNodeError> {
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
            Self::SimpleExpansion(x) => write!(f, "simple_expansion: ${0}", x.name),
            Self::VarRound(x) => write!(f, "var_round: $({x})"),
            Self::Concatenation(v) => {
                write!(f, "concatenation: ")?;
                for a in v {
                    write!(f, "{a}")?;
                }
                Ok(())
            }
            Self::String(s) => write!(f, "string: {s}"),
            Self::Expression(x) => write!(f, "expression: {x}"),
            Self::Group => write!(f, "group"),
            Self::UnderscoreIdent(x) => write!(f, "underscore_ident: {x}"),
            Self::IndexedIdent(x, index) => write!(f, "indexed_ident: {x}[{index}]"),
            Self::Word(x) => write!(f, "word: {x}"),
            Self::Error => write!(f, "ERROR"),
        }
    }
}

#[cfg(test)]
mod test {
    use pretty_assertions::assert_eq;

    use crate::utils;

    use super::*;

    #[test]
    fn simple_expansion() {
        let contents = "test_cmd $x\n";

        let tree = utils::testing::parse(contents);

        let expansion_node = tree
            .root_node()
            .child(0)
            .unwrap()
            .child(1)
            .unwrap()
            .child(0)
            .unwrap();

        dbg!(expansion_node.to_sexp());

        assert_eq!(
            Argument::from_node(&expansion_node, contents).expect("Should parse"),
            Argument {
                kind: ArgumentKind::SimpleExpansion(crate::ast::Ident {
                    name: "x".into(),
                    ident_type: crate::ast::IdentType::Variable,
                    span: ((0, 10)..(0, 11)).into(),
                }),
                span: ((0, 9)..(0, 11)).into()
            }
        )
    }

    #[test]
    fn indexed_arg() {
        let contents = "test_cmd v_x[1]\n";

        let tree = utils::testing::parse(contents);

        let expansion_node = tree
            .root_node()
            .child(0)
            .expect("command")
            .child(1)
            .expect("argument list")
            .child(0)
            .expect("indexed arg");

        dbg!(expansion_node.to_sexp());

        assert_eq!(
            Argument::from_node(&expansion_node, contents).expect("Should parse"),
            Argument {
                kind: ArgumentKind::IndexedIdent(
                    crate::ast::Ident {
                        name: "x".into(),
                        ident_type: crate::ast::IdentType::Variable,
                        span: ((0, 11)..(0, 12)).into(),
                    },
                    Index::Int(1)
                ),

                span: ((0, 9)..(0, 15)).into()
            }
        )
    }
}

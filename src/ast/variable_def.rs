use tree_sitter::Node;

use crate::{ast::ArgumentKind, ast::Ident, spans::Span, utils::tree_sitter_helpers::NodeExt};

use super::{
    from_node::{FromNode, FromNodeError},
    Argument, Word,
};

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

    fn from_node(node: &Node, text: &str) -> Result<Self, Self::Error> {
        let span: Span = node.range().into();
        let mut cursor = node.walk();

        let mut children = node.children(&mut cursor);

        let variable_keyword = children
            .next()
            // Should basically be impossible
            .ok_or(Self::Error::PartialNode(
                "missing variable keyword".to_string(),
            ))?;

        debug_assert_eq!(variable_keyword.str_text(text), "variable");

        let variable_ident = children
            .next()
            // Should basically be impossible
            .ok_or(Self::Error::PartialNode(
                "missing variable identifier".to_string(),
            ))?;

        let variable_ident = Ident::new(&variable_ident, text)?;

        let variable_kind = Word::from_node(
            &children.next().ok_or(Self::Error::PartialNode(
                "missing variable style".to_string(),
            ))?,
            text,
        )?;

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

        if variable_kind.contents != "delete" && args.is_empty() {
            return Err(Self::Error::PartialNode(
                "missing arguments in variable command".to_string(),
            ));
        }

        // TODO: Ensure that the style is valid.

        // These styles expect just a single expression.
        //  "vector" | TODO: Re-add vector after support is added for expressions in vector command
        if matches!(variable_kind.contents.as_str(), "equal" | "atom") {
            if !matches!(args[0].kind, ArgumentKind::Expression(_)) {
                return Err(Self::Error::PartialNode(format!(
                    "expected expression for variable style {}",
                    variable_kind.contents
                )));
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
mod test {

    use pretty_assertions::assert_eq;

    use super::*;
    use crate::ast;
    use crate::utils::testing::parse;

    #[test]
    fn parse_variable_def() {
        let source = "variable a equal 1.0*dt";

        let tree = parse(source);

        let root_node = tree.root_node();

        use ast::expressions::BinaryOp;
        use ast::Expression as Exp;
        let expected = VariableDef {
            variable_id: Ident {
                name: "a".into(),
                ident_type: ast::IdentType::Variable,
                span: ((0, 9)..(0, 10)).into(),
            },
            variable_style: Word {
                contents: "equal".into(),
                span: ((0, 11)..(0, 16)).into(),
            },

            args: vec![Argument {
                kind: ArgumentKind::Expression(Exp::BinaryOp(
                    Box::new(Exp::Float(1.0)),
                    BinaryOp::Multiply,
                    Box::new(Exp::ThermoKeyword(Word {
                        contents: "dt".to_owned(),
                        span: ((0, 21)..(0, 23)).into(),
                    })),
                )),
                span: ((0, 17)..(0, 23)).into(),
            }],
            span: ((0, 0)..(0, 23)).into(),
        };

        let variable_node = root_node.child(0).expect("Should find child node.");
        let parsed = VariableDef::from_node(&variable_node, source);
        assert_eq!(parsed, Ok(expected));
    }

    #[test]
    fn parse_index_variable() {
        let text = "variable file_name index step4.1.atm\n";
        // let text = include_str!("../../example_input_scripts/in.variable_index");

        let tree = parse(text);

        dbg!(tree.root_node().to_sexp());
        let root_node = tree.root_node();
        let var_def_node = root_node.child(0).unwrap(); // variable node

        let expected = VariableDef {
            variable_id: Ident {
                name: "file_name".to_string(),
                ident_type: ast::IdentType::Variable,
                span: ((0, 9)..(0, 18)).into(),
            },
            variable_style: ast::Word::new("index".to_string(), (0, 19)..(0, 24)),
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

        let ast = ast::ts_to_ast(&tree, text);
        assert!(ast.is_ok());

        if let Ok(ast) = ast {
            dbg!(&ast);
            assert_eq!(ast.commands.len(), 1);
            match &ast.commands[0] {
                crate::ast::Command::VariableDef(var) => {
                    assert_eq!(*var, expected)
                }
                cmd => panic!("Unexpected command {cmd:?}"),
            }
        }
    }

    #[test]
    fn parse_incomplete_variable_def() {
        let source = "variable a equal\n";

        let tree = parse(source);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        assert!(!root_node.is_missing());

        let expected =
            FromNodeError::PartialNode("missing arguments in variable command".to_string());

        let variable_node = root_node.child(0).expect("Should find child node.");
        let parsed = VariableDef::from_node(&variable_node, source);
        assert_eq!(parsed, Err(expected));
    }
}

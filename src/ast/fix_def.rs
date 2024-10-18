use tree_sitter::Node;

use crate::{
    ast::Ident,
    spans::Span,
    styles::FixStyle,
    utils::{into_error::IntoError, tree_sitter_helpers::NodeExt},
};

use super::{
    from_node::{FromNode, FromNodeError},
    Argument, Word,
};

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

    fn from_node(node: &Node, text: &str) -> Result<Self, Self::Error> {
        let span = node.range().into();
        let mut cursor = node.walk();

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
        let fix_style = children.next().into_err()?.str_text(text).into();

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

#[cfg(test)]
mod test {
    // Allowed within the tests
    #![allow(clippy::unwrap_used)]
    #![allow(clippy::expect_used)]

    use pretty_assertions::assert_eq;

    use super::*;
    use crate::ast;
    use crate::utils::testing::parse;

    #[test]
    fn parse_fix_no_args() {
        let source = "fix NVE all nve";

        let tree = parse(source);

        let root_node = tree.root_node();
        // Lots of tedium to parsing this...
        let fix_node = dbg!(root_node.child(0)).unwrap();

        dbg!(fix_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            FixDef::from_node(&fix_node, source).unwrap(),
            FixDef {
                fix_id: ast::Ident {
                    name: "NVE".into(),
                    ident_type: ast::IdentType::Fix,
                    span: ((0, 4)..(0, 7)).into()
                },
                group_id: ast::Word {
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
    fn parse_fix_with_args() {
        let source = "fix NVT all nvt temp 1 1.5 $(100.0*dt)
";

        let tree = parse(source);

        // let ast = ts_to_ast(&tree, source_bytes);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let fix_node = dbg!(root_node.child(0)).unwrap();

        dbg!(fix_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        use ast::ArgumentKind as AK;
        assert_eq!(
            // ast.commands[0],
            FixDef::from_node(&fix_node, source).unwrap(),
            FixDef {
                fix_id: ast::Ident {
                    name: "NVT".into(),
                    ident_type: ast::IdentType::Fix,
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
                        AK::VarRound(ast::Expression::BinaryOp(
                            ast::Expression::Float(100.0).into(),
                            ast::expressions::BinaryOp::Multiply,
                            ast::Expression::ThermoKeyword(Word {
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
}

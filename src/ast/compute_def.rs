use tree_sitter::Node;

use crate::{ast::Ident, spans::Span, styles::ComputeStyle, utils::tree_sitter_helpers::NodeExt};

use super::{
    from_node::{FromNode, FromNodeError},
    Argument, Word,
};

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
    fn from_node(node: &Node, text: &str) -> Result<Self, Self::Error> {
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
            .str_text(text)
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

#[cfg(test)]
mod test {

    use pretty_assertions::assert_eq;

    use super::*;
    use crate::ast;
    use crate::utils::testing::parse;

    #[test]
    fn parse_compute_with_args() {
        let source = "compute T_hot water temp/region hot_region";

        let tree = parse(source);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        // Lots of tedium to parsing this...
        let compute_node = dbg!(root_node.child(0)).unwrap();

        dbg!(compute_node.to_sexp());
        assert_eq!(
            ComputeDef::from_node(&compute_node, source).unwrap(),
            ComputeDef {
                compute_id: ast::Ident {
                    name: "T_hot".into(),
                    ident_type: ast::IdentType::Compute,
                    span: ((0, 8), (0, 13)).into()
                },
                group_id: Word {
                    contents: "water".into(),
                    span: ((0, 14), (0, 19)).into()
                },
                compute_style: ComputeStyle::TempRegion,
                // TODO: Change to a more generic word argument
                args: vec![ast::Argument::new(
                    ast::ArgumentKind::Word("hot_region".into()),
                    ((0, 32), (0, 42))
                ),],
                span: ((0, 0), (0, 42)).into()
            }
        );
    }

    #[test]
    fn parse_incomplete_compute_def() {
        let source = "compute a ";

        let tree = parse(source);

        let root_node = tree.root_node();
        dbg!(root_node.to_sexp());
        assert!(!root_node.is_missing());

        let ast = ast::ts_to_ast(&tree, source);

        // TODO: double check more about the syntax tree.
        assert!(ast.is_err());

        let compute_node = root_node.child(0).expect("Should find child node.");
        let parsed = ComputeDef::from_node(&compute_node, source);
        use crate::ast::from_node::FromNodeError;
        assert_eq!(
            parsed,
            Err(FromNodeError::PartialNode("missing group ID".to_string()))
        );
    }
}

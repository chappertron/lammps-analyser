use crate::ast::{self, Ast};
use crate::diagnostics::{self, Issue};
use crate::spans::{Point, Span};
use crate::styles::{ComputeStyle, FixStyle, PairStyle};
use anyhow::Result;
use once_cell::sync::Lazy;
use std::fmt::Display;
use thiserror::Error;
use tree_sitter::{Query, QueryCursor, Tree};

#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("invalid {} `{}`", style_type, name)]
pub struct InvalidStyle {
    pub start: Point,
    pub end: Point,
    pub name: String,
    pub style_type: StyleType,
}

impl Issue for InvalidStyle {
    fn diagnostic(&self) -> diagnostics::Diagnostic {
        diagnostics::Diagnostic {
            name: "invalid style".to_string(),
            severity: diagnostics::Severity::Error,
            span: Span {
                start: self.start,
                end: self.end,
            },
            message: format!("invalid {}: `{}`", self.style_type, self.name),
        }
    }
}

impl From<InvalidStyle> for lsp_types::Diagnostic {
    fn from(value: InvalidStyle) -> Self {
        lsp_types::Diagnostic::new_simple(
            lsp_types::Range {
                start: value.start.into_lsp_type(),
                end: value.end.into_lsp_type(),
            },
            value.to_string(),
        )
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum StyleType {
    Fix,
    Compute,
    Pair,
    Kspace,
    Minimize,
}

impl Display for StyleType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fix => write!(f, "fix style"),
            Self::Compute => write!(f, "compute style"),
            Self::Pair => write!(f, "pair style"),
            Self::Kspace => write!(f, "kspace style"),
            Self::Minimize => write!(f, "minimize style"),
        }
    }
}
static STYLE_QUERY: Lazy<Query> = Lazy::new(|| {
    Query::new(
        tree_sitter_lammps::language(),
        "(fix (fix_style) @definition.fix) (compute (compute_style) @definition.compute) ",
    )
    .expect("Invalid TS query")
});

// TODO: Check these by using the AST
// TODO: Also check other types of style.
/// Checks the tree for different fix and compute styles and checks if they exist or not!!!
pub fn check_styles(ast: &Ast, tree: &Tree, text: impl AsRef<[u8]>) -> Result<Vec<InvalidStyle>> {
    let text = text.as_ref();

    let mut query_cursor = QueryCursor::new();

    let matches = query_cursor.matches(&STYLE_QUERY, tree.root_node(), text);

    let mut computes_and_fixes: Vec<_> = matches
        .into_iter()
        .filter_map(|mat| {
            let style = mat.captures[0].node.utf8_text(text).ok()?;

            let style_type = match mat.captures[0].node.kind() {
                "fix_style" => StyleType::Fix,
                "compute_style" => StyleType::Compute,
                _ => unreachable!(),
            };

            // TODO: make this a match statement somehow
            if let (StyleType::Fix, FixStyle::InvalidStyle) = (&style_type, style.into()) {
                Some(InvalidStyle {
                    start: mat.captures[0].node.start_position().into(),
                    end: mat.captures[0].node.end_position().into(),
                    name: style.to_string(),
                    style_type: StyleType::Fix,
                })
            } else if let (StyleType::Compute, ComputeStyle::InvalidComputeStyle) =
                (&style_type, style.into())
            {
                Some(InvalidStyle {
                    start: mat.captures[0].node.start_position().into(),
                    end: mat.captures[0].node.end_position().into(),
                    name: style.to_string(),
                    style_type: StyleType::Compute,
                })
            } else {
                None
            }
        })
        .collect();

    // Check pair styles
    for command in &ast.commands {
        if let ast::Command::GenericCommand(cmd) = &command {
            if cmd.name.contents != "pair_style" {
                continue;
            }

            if cmd.args.is_empty() {
                // TODO: Come up with a better way of showing that the type is missing than just an
                // empty string
                // Perhaps this code is better served in the `check` code.
                computes_and_fixes.push(InvalidStyle {
                    start: cmd.end,
                    end: cmd.end,
                    name: "".to_owned(),
                    style_type: StyleType::Pair,
                });
                continue;
            }

            if let ast::ArgumentKind::Word(style) = &cmd.args[0].kind {
                if PairStyle::from(style.as_str()) == PairStyle::InvalidStyle {
                    let span = cmd.args[0].span;
                    computes_and_fixes.push(InvalidStyle {
                        start: span.start,
                        end: span.end,
                        name: style.clone(),
                        style_type: StyleType::Pair,
                    })
                }
            }
        }
    }

    Ok(computes_and_fixes)
}

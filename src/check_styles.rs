use crate::ast::{self, Ast};
use crate::diagnostics::{self, Issue};
use crate::spans::{Point, Span};
use crate::styles::{ComputeStyle, FixStyle, PairStyle};
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

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum StyleType {
    Fix,
    Compute,
    Pair,
    Kspace,
    Minimize,
}

/// Checks if fix, compute and pair styles are valid.
pub fn check_styles(ast: &Ast, tree: &Tree, text: &str) -> Vec<InvalidStyle> {
    let text = text.as_ref();

    // TODO: Check these by using the AST
    let mut invalid_styles = check_compute_and_fix_styles(tree, text);

    // Check pair styles
    // TODO: Also check other types of style.
    for command in &ast.commands {
        if let ast::Command::GenericCommand(cmd) = &command {
            if cmd.name.contents != "pair_style" {
                continue;
            }

            if cmd.args.is_empty() {
                // TODO: Come up with a better way of showing that the type is missing than just an
                // empty string
                // Perhaps this code is better served in the `check` code.
                invalid_styles.push(InvalidStyle {
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
                    invalid_styles.push(InvalidStyle {
                        start: span.start,
                        end: span.end,
                        name: style.clone(),
                        style_type: StyleType::Pair,
                    })
                }
            }
        }
    }

    invalid_styles
}

static STYLE_QUERY: Lazy<Query> = Lazy::new(|| {
    Query::new(
        &tree_sitter_lammps::LANGUAGE.into(),
        "(fix (fix_style) @definition.fix) (compute (compute_style) @definition.compute) ",
    )
    .expect("Invalid TS query")
});

fn check_compute_and_fix_styles(tree: &Tree, text: &[u8]) -> Vec<InvalidStyle> {
    let mut query_cursor = QueryCursor::new();

    let matches = query_cursor.matches(&STYLE_QUERY, tree.root_node(), text);

    matches
        .into_iter()
        .filter_map(|mat| {
            let style = mat.captures[0]
                .node
                .utf8_text(text)
                .expect("Should be valid UTF-8");

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
        .collect()
}

impl Issue for InvalidStyle {
    fn diagnostic(&self) -> diagnostics::Diagnostic {
        diagnostics::Diagnostic {
            name: "invalid style",
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

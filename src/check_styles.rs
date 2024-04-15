use crate::compute_styles::ComputeStyle;
use crate::utils::point_to_position;
use crate::{diagnostic_report::ReportSimple, fix_styles::FixStyle};
use anyhow::Result;
use owo_colors::OwoColorize;
use std::fmt::Display;
use thiserror::Error;
use tree_sitter::{Point, Query, QueryCursor, Tree};

#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("{}:{}: {} {} `{}`",start.row+1,start.column+1,"Invalid",style_type,name)]
pub struct InvalidStyle {
    pub start: Point,
    pub end: Point,
    pub name: String,
    pub style_type: StyleType,
}
impl ReportSimple for InvalidStyle {
    fn make_simple_report(&self) -> String {
        format!(
            "{}:{}: {} {} `{}`",
            self.start.row + 1,
            self.start.column + 1,
            "Invalid".bright_red(),
            self.style_type.bright_red(),
            self.name
        )
    }
}
impl From<InvalidStyle> for lsp_types::Diagnostic {
    fn from(value: InvalidStyle) -> Self {
        lsp_types::Diagnostic::new_simple(
            lsp_types::Range {
                start: point_to_position(&value.start),
                end: point_to_position(&value.end),
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

// TODO: Check these by using the AST
/// Checks the tree for different fix and compute styles and checks if they exist or not!!!
pub fn check_styles(tree: &Tree, text: impl AsRef<[u8]>) -> Result<Vec<InvalidStyle>> {
    let text = text.as_ref();
    let query = Query::new(
        tree.language(),
        "(fix (fix_style) @definition.fix) (compute (compute_style) @definition.compute) ",
    )?;

    let mut query_cursor = QueryCursor::new();

    let matches = query_cursor.matches(&query, tree.root_node(), text);

    Ok(matches
        .into_iter()
        .filter_map(|mat| {
            let style = mat.captures[0].node.utf8_text(text).ok()?;
            let style_type = match mat.captures[0].node.kind() {
                "fix_style" => StyleType::Fix,
                "compute_style" => StyleType::Compute,
                _ => unreachable!(),
            };

            if let (StyleType::Fix, FixStyle::InvalidFixStyle) = (&style_type, style.into()) {
                Some(InvalidStyle {
                    start: mat.captures[0].node.start_position(),
                    end: mat.captures[0].node.end_position(),
                    name: style.to_string(),
                    style_type: StyleType::Fix,
                })
            } else if let (StyleType::Compute, ComputeStyle::InvalidComputeStyle) =
                (&style_type, style.into())
            {
                Some(InvalidStyle {
                    start: mat.captures[0].node.start_position(),
                    end: mat.captures[0].node.end_position(),
                    name: style.to_string(),
                    style_type: StyleType::Compute,
                })
            } else {
                None
            }
        })
        .collect())
}

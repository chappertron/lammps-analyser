use std::fmt::Display;

use tree_sitter::Node;

use crate::{
    spanned_error::SpannedError,
    spans::{Point, Span},
    utils::tree_sitter_helpers::NodeExt,
};

use super::from_node::{FromNode, FromNodeError};

#[derive(Debug, Clone)]
/// Identifiers for LAMMPS fixes, computes, and variables
/// Hashing only uses the name and type, not locations
pub struct Ident {
    pub name: String,
    pub ident_type: IdentType,
    pub span: Span,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy, Default)]
pub enum IdentType {
    #[default]
    Variable,
    Fix,
    Compute,
}

impl FromNode for Ident {
    type Error = FromNodeError;

    fn from_node(node: &Node, text: &str) -> std::result::Result<Self, Self::Error> {
        let ident_type = match node.kind() {
            "fix_id" => IdentType::Fix,
            "compute_id" => IdentType::Compute,
            "variable" => IdentType::Variable,
            "f_" => IdentType::Fix,
            "c_" => IdentType::Compute,
            "v_" => IdentType::Variable,
            k => {
                return Err(FromNodeError::PartialNode(format!(
                    "invalid identifier kind {k}"
                )))
            }
        };

        let name = node.str_text(text).to_string();

        Ok(Ident {
            name,
            ident_type,
            span: node.range().into(),
        })
    }
}

impl Ident {
    /// Parses the node and text into an Ident
    /// Uses the node kind to determine the identifinder type
    /// Uses the text to extract the name of the identifier
    ///
    /// # Panics
    /// Panics if the node matches an invalid identifier type.
    pub fn new(node: &Node, text: &str) -> Result<Self, SpannedError<FromNodeError>> {
        Self::from_node(node, text).map_err(|error| SpannedError {
            error,
            span: node.range().into(),
        })
    }

    pub fn range(&self) -> Span {
        self.span
    }

    pub fn start(&self) -> Point {
        self.span.start
    }

    pub fn end(&self) -> Point {
        self.span.end
    }

    /// Display the identifier as in underscore form.
    pub fn underscore_ident(&self) -> String {
        let prefix = match self.ident_type {
            IdentType::Fix => "f_",
            IdentType::Compute => "c_",
            IdentType::Variable => "v_",
        };

        format!("{}{}", prefix, self.name)
    }
}

impl Default for Ident {
    fn default() -> Self {
        Ident {
            name: String::default(),
            ident_type: IdentType::default(),
            span: Span {
                start: Point::default(),
                end: Point::default(),
            },
        }
    }
}

impl Display for Ident {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.ident_type {
            IdentType::Fix => write!(f, "fix {}", self.name),
            IdentType::Compute => write!(f, "compute {}", self.name),
            IdentType::Variable => write!(f, "variable {}", self.name),
        }
    }
}
impl PartialEq for Ident {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name && self.ident_type == other.ident_type
    }
}

impl Eq for Ident {}

impl std::hash::Hash for Ident {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        self.ident_type.hash(state);
    }
}

impl From<&str> for IdentType {
    /// Converts from capture name to `IdentType`
    ///
    /// # Panics
    /// Panics if the string does not end with "fix", "compute", or "variable"
    fn from(value: &str) -> Self {
        if value.to_lowercase().ends_with("compute") {
            IdentType::Compute
        } else if value.to_lowercase().ends_with("fix") {
            IdentType::Fix
        } else if value.to_lowercase().ends_with("variable") {
            IdentType::Variable
        } else {
            panic!("Unknown identifier type")
        }
    }
}

impl Display for IdentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                IdentType::Fix => "fix",
                IdentType::Compute => "compute",
                IdentType::Variable => "variable",
            }
        )
    }
}

impl From<IdentType> for lsp_types::SymbolKind {
    fn from(value: IdentType) -> Self {
        match value {
            IdentType::Fix => lsp_types::SymbolKind::FUNCTION,
            IdentType::Compute => lsp_types::SymbolKind::FUNCTION,
            IdentType::Variable => lsp_types::SymbolKind::VARIABLE,
        }
    }
}

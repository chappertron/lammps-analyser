use crate::{
    ast::from_node::{FromNode, FromNodeError},
    diagnostics::Issue,
    spans::{Point, Span},
};
use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
};
use thiserror::Error;
use tree_sitter::{Node, Query, QueryCursor, Tree};

use once_cell::sync::Lazy;

pub type IdentMap = HashMap<NameAndType, SymbolDefsAndRefs>;

/// Find and store Identifiers in the `tree-sitter` Tree. Stores a `QueryCursor` for re-use
/// Symbols can be accessed through the `symbols` method.
pub struct IdentiFinder {
    cursor: QueryCursor,
    symbols: HashMap<NameAndType, SymbolDefsAndRefs>,
}

impl Debug for IdentiFinder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IdentiFinder")
            .field("symbols", &self.symbols)
            .finish()
    }
}

#[derive(Debug, Clone, Default)]
pub struct SymbolDefsAndRefs {
    /// Definition or definitions of the symbol
    defs: SymbolDef,
    /// References to the symbol
    refs: Vec<Ident>,
}

impl SymbolDefsAndRefs {
    pub fn has_def(&self) -> bool {
        !self.defs.is_none()
    }

    pub fn defs(&self) -> &SymbolDef {
        &self.defs
    }
    pub fn refs(&self) -> &[Ident] {
        &self.refs
    }
}

#[derive(Debug, Clone, Default)]
pub struct SymbolDef {
    pub defs: Vec<Ident>,
}

impl SymbolDef {
    /// Returns `true` if the defs is [`None`].
    ///
    /// [`None`]: Defs::None
    #[must_use]
    fn is_none(&self) -> bool {
        self.defs.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Ident> {
        self.defs.iter()
    }
}

static QUERY_DEF: Lazy<Query> = Lazy::new(|| {
    Query::new(
        tree_sitter_lammps::language(),
        "(fix (fix_id ) @definition.fix) 
                (compute (compute_id) @definition.compute) 
                (variable_def (variable) @definition.variable )",
    )
    .expect("Invalid query for LAMMPS TS Grammar")
});

static QUERY_REF: Lazy<Query> = Lazy::new(|| {
    Query::new(
        tree_sitter_lammps::language(),
        " (fix_id) @reference.fix  (compute_id) @reference.compute (variable) @reference.variable",
    )
    .expect("Invalid query for LAMMPS TS Grammar")
});

impl IdentiFinder {
    /// Creates a new `IdentiFinder` without searching for the symbols, leaving an empty `symbols`
    /// map
    ///
    /// Fails if the Query used internally is invalid.
    pub fn new_no_parse() -> Self {
        IdentiFinder {
            cursor: QueryCursor::new(),
            symbols: HashMap::new(),
        }
    }

    /// Create a new `Identifinder` and search for symbols.
    pub fn new(tree: &Tree, text: impl AsRef<[u8]>) -> Result<Self, FromNodeError> {
        let text = text.as_ref();
        let mut i = Self::new_no_parse();
        i.find_symbols(tree, text)?;
        Ok(i)
    }

    pub fn find_symbols(
        &mut self,
        tree: &Tree,
        text: &[u8],
    ) -> Result<&HashMap<NameAndType, SymbolDefsAndRefs>, FromNodeError> {
        let captures = self.cursor.captures(&QUERY_DEF, tree.root_node(), text);

        // TODO: Clear the symbols in a smarter way, perhaps only removing old ones
        // Do this for an incremental method
        self.symbols.clear();

        for (mtch, _cap_id) in captures {
            let node = mtch.captures[0].node;

            let ident = Ident::new(&node, text)?;

            let name_and_type = NameAndType {
                name: ident.name.clone(),
                ident_type: ident.ident_type,
            };

            let entry = self.symbols.entry(name_and_type).or_default();

            entry.defs.defs.push(ident);
        }
        // references
        let captures = self.cursor.captures(&QUERY_REF, tree.root_node(), text);

        for (mtch, _cap_id) in captures {
            let node = mtch.captures[0].node;

            let ident = Ident::new(&node, text)?;

            let name_and_type = NameAndType {
                name: ident.name.clone(),
                ident_type: ident.ident_type,
            };

            let entry = self.symbols.entry(name_and_type).or_default();
            entry.refs.push(ident);
        }

        Ok(&self.symbols)
    }

    /// Check for symbols that have been used without being defined
    /// TODO: Does not currently verify order or if the fix has been deleted at point of definition
    /// TODO: Return the Vec or an option, not an error
    pub fn check_symbols(&self) -> Result<(), Vec<UndefinedIdent>> {
        let undefined_fixes: Vec<_> = self
            .symbols
            .iter()
            .filter_map(|(_k, v)| {
                // TODO: Double check this? Seems opposite to what
                if v.defs.is_none() {
                    Some(v.refs.iter())
                } else {
                    None
                }
            })
            .flatten()
            .map(|x| UndefinedIdent { ident: x.clone() })
            .collect();

        if undefined_fixes.is_empty() {
            Ok(())
        } else {
            Err(undefined_fixes)
        }
    }

    pub fn symbols(&self) -> &HashMap<NameAndType, SymbolDefsAndRefs> {
        &self.symbols
    }
}
/// Finds all the unused variables in the input file.
/// TODO: Convert to a warning type
pub fn unused_variables(map: &HashMap<NameAndType, SymbolDefsAndRefs>) -> Vec<UnusedIdent> {
    map.iter()
        .filter_map(|(k, v)| {
            // TODO: don't include definitions as references???
            if v.refs().len() == v.defs().defs.len() && k.ident_type == IdentType::Variable {
                Some(v.refs())
            } else {
                None
            }
        })
        .flatten()
        // TODO: BOO CLONE
        .map(|x| UnusedIdent { ident: x.clone() })
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("unused {} `{}`", ident.ident_type, ident.name)]
pub struct UnusedIdent {
    pub ident: Ident,
}

impl From<UnusedIdent> for lsp_types::Diagnostic {
    fn from(value: UnusedIdent) -> Self {
        lsp_types::Diagnostic {
            message: format!("unused {} `{}`", value.ident.ident_type, value.ident.name),
            range: value.ident.range().into_lsp_types(),
            severity: Some(lsp_types::DiagnosticSeverity::WARNING),
            ..Default::default()
        }
    }
}

impl From<Ident> for UnusedIdent {
    fn from(ident: Ident) -> Self {
        Self { ident }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct NameAndType {
    pub name: String,
    pub ident_type: IdentType,
}

#[derive(Debug, Clone)]
/// Identifiers for LAMMPS fixes, computes, and variables
/// Hashing only uses the name and type, not locations
pub struct Ident {
    pub name: String,
    pub ident_type: IdentType,
    pub span: Span,
}

impl FromNode for Ident {
    type Error = FromNodeError;

    fn from_node(node: &Node, text: impl AsRef<[u8]>) -> std::result::Result<Self, Self::Error> {
        let text = text.as_ref();
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

        let name = node.utf8_text(text)?.to_string();

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
    pub fn new(node: &Node, text: &[u8]) -> Result<Self, FromNodeError> {
        Self::from_node(node, text)
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

impl Hash for Ident {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        self.ident_type.hash(state);
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy, Default)]
pub enum IdentType {
    #[default]
    Variable,
    Fix,
    Compute,
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

#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("undefined {} `{}`", ident.ident_type,ident.name)]
pub struct UndefinedIdent {
    pub ident: Ident,
}
impl Issue for UnusedIdent {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        let name = "unused identifier";
        crate::diagnostics::Diagnostic {
            name: name.to_string(),
            severity: crate::diagnostics::Severity::Warning,
            span: self.ident.span,
            message: self.to_string(),
        }
    }
}
impl Issue for UndefinedIdent {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        let name = "undefined identifier";
        crate::diagnostics::Diagnostic {
            name: name.to_string(),
            severity: crate::diagnostics::Severity::Error,
            span: self.ident.span,
            message: self.to_string(),
        }
    }
}

impl From<UndefinedIdent> for lsp_types::Diagnostic {
    fn from(value: UndefinedIdent) -> Self {
        lsp_types::Diagnostic::new_simple(value.ident.range().into_lsp_types(), value.to_string())
    }
}

#[cfg(test)]
mod tests {}

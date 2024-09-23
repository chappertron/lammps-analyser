//! Definies the `Identifinder` type which finds and validates definitions and references.

use crate::ast::{Ident, IdentType};
use crate::spanned_error::SpannedError;
use crate::{ast::from_node::FromNodeError, diagnostics::Issue};
use std::{collections::HashMap, fmt::Debug};
use thiserror::Error;
use tree_sitter::{Query, QueryCursor, Tree};

use once_cell::sync::Lazy;

pub type IdentMap = HashMap<NameAndType, SymbolDefsAndRefs>;

/// Find and store Identifiers in the `tree-sitter` Tree. Stores a `tree_sitter::QueryCursor` for re-use
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
    /// Returns the `false` if the list of identifiers are empty
    #[must_use]
    fn is_none(&self) -> bool {
        self.defs.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Ident> {
        self.defs.iter()
    }
}

/// Query for identifier definitions
static QUERY_DEF: Lazy<Query> = Lazy::new(|| {
    Query::new(
        tree_sitter_lammps::language(),
        "(fix (fix_id ) @definition.fix) 
                (compute (compute_id) @definition.compute) 
                (variable_def (variable) @definition.variable )",
    )
    .expect("Invalid query for LAMMPS TS Grammar")
});

/// Query for identifier references
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
    pub fn new(tree: &Tree, text: &str) -> Result<Self, FromNodeError> {
        let text = text.as_ref();
        let mut i = Self::new_no_parse();
        i.find_symbols(tree, text)?;
        Ok(i)
    }

    pub fn find_symbols(
        &mut self,
        tree: &Tree,
        text: &str,
    ) -> Result<&HashMap<NameAndType, SymbolDefsAndRefs>, SpannedError<FromNodeError>> {
        let captures = self
            .cursor
            .captures(&QUERY_DEF, tree.root_node(), text.as_bytes());

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
        let captures = self
            .cursor
            .captures(&QUERY_REF, tree.root_node(), text.as_bytes());

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
    ///
    /// If invalid references are found, an error is returned.
    pub fn check_symbols(&self) -> Result<(), Vec<UndefinedIdent>> {
        // TODO: Does not currently verify order or if the fix has been deleted at point of definition

        let undefined_fixes: Vec<_> = self
            .symbols
            .iter()
            .filter_map(|(_k, v)| {
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
///
/// Note: there is no equivalent for fixes and computes as these generally have sideffects
/// TODO: add an equivalent for computes. Current blocker is false negatives.
pub fn unused_references(map: &HashMap<NameAndType, SymbolDefsAndRefs>) -> Vec<UnusedIdent> {
    map.iter()
        .filter_map(|(k, v)| {
            // TODO: don't include definitions as references
            if v.refs().len() == v.defs().defs.len() && matches!(k.ident_type, IdentType::Variable)
            // | IdentType::Compute
            {
                Some(v.refs())
            } else {
                None
            }
        })
        .flatten()
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

#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("undefined {} `{}`", ident.ident_type,ident.name)]
pub struct UndefinedIdent {
    pub ident: Ident,
}
impl Issue for UnusedIdent {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        let name = "unused identifier";
        crate::diagnostics::Diagnostic {
            name,
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
            name,
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

use anyhow::Result;
use owo_colors::OwoColorize;
use std::{
    collections::HashSet,
    fmt::{Debug, Display},
    hash::Hash,
};
use thiserror::Error;
use tree_sitter::{Node, Point, Query, QueryCursor, Tree};

use crate::{diagnostic_report::ReportSimple, utils::point_to_position};

pub struct IdentiFinder {
    pub query_def: Query,
    pub query_ref: Query,
    cursor: QueryCursor,
    ident_defs: HashSet<Ident>,
    ident_refs: Vec<Ident>,
}

impl IdentiFinder {
    pub fn new(tree: &Tree, text: &[u8]) -> Result<Self> {
        let query_def = Query::new(
            tree_sitter_lammps::language(),
            " (fix (fix_id ) @definition.fix) (compute (compute_id) @definition.compute) (variable_def (variable) @definition.variable )",
        )?;

        let query_ref = Query::new(
            tree_sitter_lammps::language(),
            " (fix_id) @reference.fix  (compute_id) @reference.compute (variable) @reference.variable",
        )?;

        let mut i = IdentiFinder {
            query_def,
            query_ref,
            cursor: QueryCursor::new(),
            ident_defs: HashSet::new(),
            ident_refs: Vec::new(),
        };

        i.find_refs(tree, text)?;
        i.find_defs(tree, text)?;

        Ok(i)
    }

    pub fn find_refs(&mut self, tree: &Tree, text: &[u8]) -> Result<&Vec<Ident>> {
        let captures = self
            .cursor
            .captures(&self.query_ref, tree.root_node(), text);

        self.ident_refs.extend(captures.map(|(mtch, _cap_id)| {
            let node = mtch.captures[0].node;

            // TODO Properly handle errors
            Ident::new(&node, text).unwrap()
        }));

        Ok(&self.ident_refs)
    }

    pub fn find_defs(&mut self, tree: &Tree, text: &[u8]) -> Result<&HashSet<Ident>> {
        let captures = self
            .cursor
            .captures(&self.query_def, tree.root_node(), text);

        self.ident_defs.extend(captures.map(|(mtch, _cap_id)| {
            let node = mtch.captures[0].node;

            // TODO perhaps better to check for this in another way?
            Ident::new(&node, text).unwrap()
        }));

        Ok(&self.ident_defs)
    }

    /// Check for symbols that have been used without being defined
    /// TODO: Does not currently verify order or if the fix has been deleted at point of definition
    /// TODO Return the Vec or an option, not an error
    pub fn check_symbols(&mut self) -> Result<(), Vec<UndefinedIdent>> {
        let used_fixes: HashSet<_> = self.ident_refs.iter().collect();

        let undefined_fixes: Vec<_> = used_fixes
            .difference(&self.ident_defs.iter().collect())
            .map(|&x| UndefinedIdent { ident: x.clone() })
            .collect();
        assert_ne!(undefined_fixes.len(), used_fixes.len());
        if undefined_fixes.is_empty() {
            Ok(())
        } else {
            Err(undefined_fixes)
        }
    }
}

#[derive(Debug, Clone)]
/// Identifiers for LAMMPS fixes, computes, and variables
/// Hashing only uses the name and type, not locations
pub struct Ident {
    pub name: String,
    pub ident_type: IdentType,
    pub start: Point,
    pub end: Point,
    pub start_byte: usize,
    pub end_byte: usize,
}
impl Ident {
    /// Parses the node and text into an Ident
    /// Uses the node kind to determine the identifinder type
    /// Uses the text to extract the name of the identifier
    pub fn new(node: &Node, text: &[u8]) -> Result<Self> {
        let name = node.utf8_text(text)?.to_string();

        let ident_type = match node.kind() {
            "fix_id" => IdentType::Fix,
            "compute_id" => IdentType::Compute,
            "variable" => IdentType::Variable,
            x => panic!("Unknown identifier type {x}"), // TODO Make this not panic
        };

        Ok(Ident {
            name,
            ident_type,
            start: node.start_position(),
            end: node.end_position(),
            start_byte: node.start_byte(),
            end_byte: node.end_byte(),
        })
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

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub enum IdentType {
    Fix,
    Variable,
    Compute,
}

impl From<&str> for IdentType {
    /// Converts from capture name to IdentType
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

#[derive(Debug, Error, Clone)]
#[error("{}:{}: {} {} `{}`",ident.start.row+1,ident.start.column+1,"Undefined",ident.ident_type,ident.name)]
pub struct UndefinedIdent {
    pub ident: Ident,
}

impl ReportSimple for UndefinedIdent {
    fn make_simple_report(&self) -> String {
        format!(
            "{}:{}: {} {} `{}`",
            self.ident.start.row + 1,
            self.ident.start.column + 1,
            "Undefined".bright_red(),
            self.ident.ident_type.bright_red(),
            self.ident.name
        )
    }
}

impl From<UndefinedIdent> for lsp_types::Diagnostic {
    fn from(value: UndefinedIdent) -> Self {
        lsp_types::Diagnostic::new_simple(
            lsp_types::Range {
                start: point_to_position(&value.ident.start),
                end: point_to_position(&value.ident.end),
            },
            value.to_string(),
        )
    }
}

// use crate::diagnostic_report::ReportDiagnostic;
// use ariadne::{Label, ReportKind,Span};

// impl<S:Span> ReportDiagnostic<S> for UndefinedIdent {
//     fn make_report(&self, source_id: &str) -> ariadne::Report<'_,S>
//     {
//         ariadne::Report::build(ReportKind::Error, source_id, self.ident.start_byte)
//             .with_label(
//                 Label::new((source_id, self.ident.start_byte..self.ident.end_byte))
//                     .with_message(format!("{}", self)),
//             )
//             .finish()
//     }
// }

#[cfg(test)]
mod tests {}

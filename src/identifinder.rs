use anyhow::Result;
use owo_colors::OwoColorize;
use std::{
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
};
use thiserror::Error;
use tree_sitter::{Node, Point, Query, QueryCursor, Range, Tree};

use crate::{diagnostic_report::ReportSimple, utils::point_to_position};

pub struct IdentiFinder {
    pub query_def: Query,
    pub query_ref: Query,
    cursor: QueryCursor,
    // ident_defs: Vec<SymbolDef>,
    // ident_refs: Vec<Ident>,
    symbols: HashMap<NameAndType, SymbolDefsAndRefs>,
}

impl Debug for IdentiFinder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IdentiFinder")
            .field("query_def", &self.query_def)
            .field("query_ref", &self.query_ref)
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

// #[derive(Debug, Clone, Default)]
// pub enum SymbolDef {
//     /// Single Definition
//     Single(Ident),
//     /// Multiple Defintions
//     Multiple(Vec<Ident>),
//     /// No Definition
//     #[default]
//     None,
// }

#[derive(Debug, Clone, Default)]
pub struct SymbolDef {
    defs: Vec<Ident>,
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

// TODO make this into a new data structure, holding definitions and all references in one hashmap
// Use `DashMap` for easier use with the LSP Server???

impl IdentiFinder {
    pub fn new_no_parse() -> Result<Self> {
        let query_def = Query::new(
            tree_sitter_lammps::language(),
            "(fix (fix_id ) @definition.fix) 
                (compute (compute_id) @definition.compute) 
                (variable_def (variable) @definition.variable )",
        )?;

        let query_ref = Query::new(
            tree_sitter_lammps::language(),
            " (fix_id) @reference.fix  (compute_id) @reference.compute (variable) @reference.variable",
        )?;

        let i = IdentiFinder {
            query_def,
            query_ref,
            cursor: QueryCursor::new(),
            // ident_defs: HashSet::new(),
            // ident_refs: Vec::new(),
            symbols: HashMap::new(),
        };

        // i.find_defs(tree, text)?;
        // i.find_refs(tree, text)?;

        Ok(i)
    }

    pub fn new(tree: &Tree, text: &[u8]) -> Result<Self> {
        let mut i = Self::new_no_parse()?;
        i.find_symbols(tree, text)?;
        Ok(i)
    }

    pub fn find_symbols(
        &mut self,
        tree: &Tree,
        text: &[u8],
    ) -> Result<&HashMap<NameAndType, SymbolDefsAndRefs>> {
        let captures = self
            .cursor
            .captures(&self.query_def, tree.root_node(), text);

        // TODO Clear the symbols in a smarter way, perhaps only removing old ones
        // Do this for an incremental method
        self.symbols.clear();

        for (mtch, _cap_id) in captures {
            let node = mtch.captures[0].node;

            let ident = Ident::new(&node, text)?;

            let name_and_type = NameAndType {
                name: ident.name.clone(),
                ident_type: ident.ident_type.clone(),
            };

            let entry = self.symbols.entry(name_and_type).or_default();
            // match &mut entry.defs {
            //     SymbolDef::Multiple(v) => {
            //         v.push(ident);
            //     }
            //     SymbolDef::Single(ident_0) => {
            //         // `std::mem::take` done to avoid clone. Might not be better because a new default is
            //         // still created???
            //         entry.defs = SymbolDef::Multiple(vec![std::mem::take(ident_0), ident]);
            //     }
            //     SymbolDef::None => {
            //         entry.defs = SymbolDef::Single(ident);
            //     }
            // }
            entry.defs.defs.push(ident);
        }
        // references
        let captures = self
            .cursor
            .captures(&self.query_ref, tree.root_node(), text);

        for (mtch, _cap_id) in captures {
            let node = mtch.captures[0].node;

            let ident = Ident::new(&node, text)?;

            let name_and_type = NameAndType {
                name: ident.name.clone(),
                ident_type: ident.ident_type.clone(),
            };

            let entry = self.symbols.entry(name_and_type).or_default();
            entry.refs.push(ident);
        }

        Ok(&self.symbols)
    }

    // pub fn find_refs(&mut self, tree: &Tree, text: &[u8]) -> Result<&Vec<Ident>> {
    //     let captures = self
    //         .cursor
    //         .captures(&self.query_ref, tree.root_node(), text);

    //     self.ident_refs.extend(captures.map(|(mtch, _cap_id)| {
    //         let node = mtch.captures[0].node;

    //         // TODO Properly handle errors
    //         Ident::new(&node, text).unwrap()
    //     }));

    //     Ok(&self.ident_refs)
    // }

    // pub fn find_defs(&mut self, tree: &Tree, text: &[u8]) -> Result<&HashSet<Ident>> {
    //     let captures = self
    //         .cursor
    //         .captures(&self.query_def, tree.root_node(), text);

    //     self.ident_defs.extend(captures.map(|(mtch, _cap_id)| {
    //         let node = mtch.captures[0].node;

    //         // TODO perhaps better to check for this in another way?
    //         Ident::new(&node, text).unwrap()
    //     }));

    //     Ok(&self.ident_defs)
    // }

    /// Check for symbols that have been used without being defined
    /// TODO: Does not currently verify order or if the fix has been deleted at point of definition
    /// TODO Return the Vec or an option, not an error
    pub fn check_symbols(&self) -> Result<(), Vec<UndefinedIdent>> {
        let undefined_fixes: Vec<_> = self
            .symbols
            .iter()
            .filter_map(|(_k, v)| {
                if !v.defs.is_none() {
                    None
                } else {
                    Some(v.refs.iter())
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
        // let used_fixes: HashSet<_> = self.ident_refs.iter().collect();

        // let undefined_fixes: Vec<_> = used_fixes
        //     .difference(&self.ident_defs.iter().collect())
        //     .map(|&x| UndefinedIdent { ident: x.clone() })
        //     .collect();
        // assert_ne!(undefined_fixes.len(), used_fixes.len());
        // if undefined_fixes.is_empty() {
        //     Ok(())
        // } else {
        //     Err(undefined_fixes)
        // }
        //
        // todo!("New version of the checking symbols!!")
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
            // TODO don't include defintions as references???
            if v.refs().len() == v.defs().defs.len() && k.ident_type == IdentType::Variable {
                Some(v.refs())
            } else {
                None
            }
        })
        .flatten()
        // BOO CLONE
        .map(|x| UnusedIdent { ident: x.clone() })
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error(
            "{}:{}: {} {} `{}`",
            ident.start.row + 1,
            ident.start.column + 1,
            "Unused",
            ident.ident_type,
            ident.name
)]
pub struct UnusedIdent {
    pub ident: Ident,
}

impl ReportSimple for UnusedIdent {
    fn make_simple_report(&self) -> String {
        format!(
            "{}:{}: {} {} `{}`",
            self.ident.start.row + 1,
            self.ident.start.column + 1,
            "Unused".bright_yellow(),
            self.ident.ident_type.bright_yellow(),
            self.ident.name
        )
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

#[derive(Debug, Clone, Default)]
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

    pub fn range(&self) -> Range {
        Range {
            start_byte: self.start_byte,
            end_byte: self.end_byte,
            start_point: self.start,
            end_point: self.end,
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

impl From<IdentType> for lsp_types::SymbolKind {
    fn from(value: IdentType) -> Self {
        match value {
            IdentType::Fix => lsp_types::SymbolKind::FUNCTION,
            IdentType::Compute => lsp_types::SymbolKind::FUNCTION,
            IdentType::Variable => lsp_types::SymbolKind::VARIABLE,
        }
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

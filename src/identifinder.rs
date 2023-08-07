use anyhow::Result;
use std::{collections::HashSet, hash::Hash};
use tree_sitter::{Node, Point, Query, Tree, QueryCursor};

pub struct IdentiFinder {
    pub query_def: Query,
    pub query_ref: Query,
    cursor : QueryCursor,
    ident_defs: HashSet<Ident>,
    ident_refs: Vec<Ident>,
}

impl IdentiFinder {
    pub fn new() -> Result<Self> {
        let query_def = Query::new(
            tree_sitter_lammps::language(),
            " (fix (fix_id ) @definition.fix) (compute (compute_id) @definition.compute)",
        )?;

        let query_ref = Query::new(
            tree_sitter_lammps::language(),
            " (fix_id) @reference.fix  (compute_id) @reference.compute",
        )?;

        Ok( IdentiFinder {  query_def, query_ref,cursor:QueryCursor::new(),ident_defs: HashSet::new(), ident_refs: Vec::new(),})
    }

    pub fn find_refs(&mut self, tree: &Tree, text: &[u8]) -> Result<&Vec<Ident>> {
       
        let captures = self.cursor.captures(&self.query_ref, tree.root_node(), text);
        
        self.ident_refs.extend(captures.map(|(mtch,cap_id)| {
            
            
            todo!("find_refs");
            todo!("Convert capture name into IdentType");

        } ));

        Ok(&self.ident_refs)
    }


    pub fn find_defs(&mut self, tree: &Tree, text: &[u8]) -> Result<&HashSet<Ident>> {
       
        let captures = self.cursor.captures(&self.query_ref, tree.root_node(), text);

        todo!("find_defs");
        Ok(&self.ident_defs)
    }


}

#[derive(Debug)]
/// Identifiers for LAMMPS fixes, computes, and variables
/// Hashing only uses the name and type, not locations
pub struct Ident {
    pub name: String,
    pub ident_type: IdentType,
    pub start: Point,
    pub end: Point,
}
impl Ident {
    pub fn new(node: &Node, text: &[u8]) -> Result<Self> {
        let name = node.utf8_text(text)?.to_string();

        let ident_type = match node.kind() {
            "fix_id" => IdentType::Fix,
            "compute_id" => IdentType::Compute,
            "var" => IdentType::Variable,
            _ => panic!("Unknown identifier type"), // TODO Make this not panic
        };

        Ok(Ident {
            name,
            ident_type,
            start: node.start_position(),
            end: node.end_position(),
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

#[derive(Debug, PartialEq, Eq, Hash)]
pub enum IdentType {
    Fix,
    Variable,
    Compute,
}

impl From<&str> for IdentType {
   fn from(value: &str) -> Self {
        todo!("Match these!!!!") 
   } 
}


#[cfg(test)]
mod tests {
    use super::IdentiFinder;


    #[test]
    fn new() {
       let identifinder = IdentiFinder::new().unwrap();
        assert!(false)
    }

}

use tree_sitter::Tree;

use crate::{
    ast::{Ast, Command, NamedCommand},
    identifinder::IdentMap,
};

/// Lint against fixes being redefined or deleted before they are run
/// Finds all fixes that meet this criteria.
pub fn fix_redef_before_run(tree: &Tree, text: &[u8], ast: &Ast, indents: &IdentMap) {
    // Find fixes with multiple definition

    for (k, v) in indents {
        if v.defs().defs.len() > 1 {
            dbg!(k);
            dbg!(v.defs());
        }
    }

    for command in &ast.commands {
        dbg!(command);
        match command {
            Command::NamedCommand(NamedCommand::Run) => {
                dbg!("run");
            }
            _ => (),
        }
    }
}

#[cfg(test)]
mod tests {
    use tree_sitter::Parser;

    use crate::{ast::ts_to_ast, identifinder::IdentiFinder};

    use super::*;
    fn setup_parser() -> Parser {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        parser
    }
    #[test]
    fn test_multiple_fix_defs() {
        let mut parser = setup_parser();
        let text = b"fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 10000";

        let tree = parser.parse(text, None).unwrap();

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents.find_symbols(&tree, text);

        fix_redef_before_run(&tree, text, &ast, idents.symbols());

        unimplemented!()
    }
}

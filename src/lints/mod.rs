use crate::ast::FixDef;
use crate::{
    ast::{Ast, CommandType, NamedCommand},
    identifinder::{Ident, IdentMap},
};

/// If a fix is redefined before it is run, the first definition is useless.
#[derive(Debug)]
pub struct MultiplyDefinedBeforeRun<'a> {
    pub first_def: &'a Ident,
    pub other_defs: Vec<&'a Ident>,
}

/// Lint against fixes being redefined or deleted before they are run
/// Finds all fixes that meet this criteria.
pub fn fix_redef_before_run<'a>(
    ast: &'a Ast,
    indents: &'a IdentMap,
) -> Result<(), Vec<MultiplyDefinedBeforeRun<'a>>> {
    // Find fixes with multiple definition

    // for (k, v) in indents {
    //     if v.defs().defs.len() > 1 {
    //         dbg!(k);
    //         dbg!(v.defs());
    //     }
    // }

    let multiply_defined = indents.values().filter(|v| v.defs().defs.len() > 1);

    // TODO: Exit early if no multiply defined commands found
    // could have negative performance impact, however, if collecting vec is slow.

    // // Find run nodes
    // let run_commands = ast.commands.iter().filter(|command| {
    //     matches!(
    //         command.command_type,
    //         CommandType::NamedCommand(NamedCommand::Run)
    //     )
    // });

    // Count the number of times a fix is defined before the next run command
    let mut redefined = vec![];

    for v in multiply_defined {
        let mut defs = v.defs().defs.iter();
        let first_def = defs.next().unwrap();

        // Skip until we find first_def
        let commands_iter = ast
            .commands
            .iter()
            // Check if definition is in range, because definiton uses identifier location.
            .skip_while(|cmd| !(first_def.range()).in_range(&cmd.range()))
            .skip(1);

        // Append to bad commands list if second def found

        // Take until the next run command
        let other_defs: Vec<_> =
            commands_iter
                .take_while(|cmd| {
                    !matches!(
                        cmd.command_type,
                        CommandType::NamedCommand(NamedCommand::Run)
                    )
                })
                .filter_map(|cmd| match &cmd.command_type {
                    // TODO: Add compute support too
                    CommandType::NamedCommand(NamedCommand::Fix(FixDef {
                        fix_id: ident, ..
                    })) if ident == first_def => Some(ident),
                    _ => None,
                })
                .collect();

        if !other_defs.is_empty() {
            redefined.push(MultiplyDefinedBeforeRun {
                first_def,
                other_defs,
            });
        }

        // let next_nearest_run_command  =
    }

    if redefined.is_empty() {
        Ok(())
    } else {
        Err(redefined)
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
                             run 10000
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = parser.parse(text, None).unwrap();

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents
            .find_symbols(&tree, text)
            .expect("Failed to Find Symbols");

        if let Err(v) = fix_redef_before_run(&ast, idents.symbols()) {
            assert_eq!(v.len(), 1);
            dbg!(v);
        } else {
            panic!("Error was expected")
        }
    }

    #[test]
    fn redefined_fix_after_run() {
        let mut parser = setup_parser();
        let text = b"fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 10000
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = parser.parse(text, None).unwrap();

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents
            .find_symbols(&tree, text)
            .expect("Failed to Find Symbols");

        assert!(fix_redef_before_run(&ast, idents.symbols()).is_ok());
    }

    #[test]
    fn redefined_twice_after_first_run() {
        let mut parser = setup_parser();
        let text = b"fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 10000
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = parser.parse(text, None).unwrap();

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents
            .find_symbols(&tree, text)
            .expect("Failed to Find Symbols");

        if let Err(v) = fix_redef_before_run(&ast, idents.symbols()) {
            assert_eq!(v.len(), 1);
            dbg!(v);
        } else {
            panic!("Error was expected")
        }
    }
}

use std::collections::HashMap;

use crate::diagnostics::Issue;
use crate::spans::Span;
use crate::{
    ast::Ast,
    identifinder::{Ident, IdentMap},
};

/// If a fix is redefined before it is run, the first definition is useless.
#[derive(Debug, PartialEq)]
struct MultiplyDefinedBeforeRun<'a> {
    defs: Vec<&'a Ident>,
}

#[derive(Debug)]
pub struct RedfinedIdent<'a>(&'a Ident);

pub fn redefined_identifiers<'a>(
    ast: &'a Ast,
    idents: &'a IdentMap,
) -> impl Iterator<Item = RedfinedIdent<'a>> {
    fix_redef_before_run(ast, idents)
        .into_iter()
        // Don't show the warning on the first ident.
        // TODO: show an Information for it though
        .flat_map(|x| x.defs.into_iter().skip(1).map(RedfinedIdent))
}

/// Lint against fixes being redefined or deleted before they are run
/// Finds all fixes that meet this criteria.
fn fix_redef_before_run<'a>(
    ast: &'a Ast,
    indents: &'a IdentMap,
) -> Vec<MultiplyDefinedBeforeRun<'a>> {
    // All the commands that have multiple definitions.
    let multiply_defined = indents.values().filter(|v| v.defs().defs.len() > 1);

    let run_blocks: Vec<_> = ast.find_run_blocks().collect();

    // Count the number of times a fix is defined before the next run command
    let mut redefined = vec![];

    for v in multiply_defined {
        let defs = v.defs();

        let mut defs_by_blocks: HashMap<Span, Vec<&Ident>> = HashMap::new();

        for block in &run_blocks {
            for def in defs.iter() {
                if block.contains_span(&def.range()) {
                    defs_by_blocks.entry(*block).or_default().push(def);
                }
            }
        }

        for (_span, defs) in defs_by_blocks {
            if defs.len() > 1 {
                redefined.push(MultiplyDefinedBeforeRun { defs })
            }
        }
    }

    redefined
}

impl Issue for RedfinedIdent<'_> {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        let ident = self.0;
        crate::diagnostics::Diagnostic {
            name: "redefined before `run` command".to_string(),
            severity: crate::diagnostics::Severity::Warning,
            span: ident.span,
            message: format![
                "{} `{}` defined multiple times before run command",
                ident.ident_type, ident.name
            ],
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use crate::{ast::ts_to_ast, identifinder::IdentiFinder, utils};

    use super::*;

    #[test]
    fn test_multiple_fix_defs() {
        let text = "fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 10000
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = utils::testing::parse(text);

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents
            .find_symbols(&tree, text)
            .expect("Failed to Find Symbols");

        let v = fix_redef_before_run(&ast, idents.symbols());
        assert_eq!(v.len(), 1);
        dbg!(v);
    }

    #[test]
    fn redefined_fix_after_run() {
        let text = "fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 10000
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)";

        let tree = utils::testing::parse(text);

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents
            .find_symbols(&tree, text)
            .expect("Failed to Find Symbols");

        assert_eq!(fix_redef_before_run(&ast, idents.symbols()), vec![]);
    }

    #[test]
    // #[ignore = "Lint being tested is incomplete"]
    // TODO: Finish the lint so the test passes
    fn redefined_twice_after_first_run() {
        let text = "fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 10000
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             fix NVT all nvt temp 1 1.5 $(100.0*dt)
                             run 20000
";

        let tree = utils::testing::parse(text);

        let ast = ts_to_ast(&tree, text).unwrap();

        let mut idents = IdentiFinder::new(&tree, text).unwrap();
        idents
            .find_symbols(&tree, text)
            .expect("Failed to Find Symbols");

        let v = fix_redef_before_run(&ast, idents.symbols());
        assert_eq!(v.len(), 1);
        dbg!(v);
    }
}

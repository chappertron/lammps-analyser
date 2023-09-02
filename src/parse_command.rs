// Parse LAMMPS commands into a struct using Clap?
// TODO I am not sure this will work. Primarily because there
// are no -- flags. I think I will need to use a custom parser :(

// use clap::Parser;

use thiserror::Error;

use crate::fix_styles::FixStyle;

#[derive(Debug, Error)]
pub enum InvalidArguments {
    #[error("Incorrect number of arguments: {provided} expected {expected}")]
    IncorrectNumberArguments { provided: usize, expected: usize },
}

/// Parse
pub fn parse_fix(args: &[&str]) -> Result<(), InvalidArguments> {
    if args.len() < 3 {
        return Err(InvalidArguments::IncorrectNumberArguments {
            provided: args.len(),
            expected: 3,
        });
    }

    // TODO use these
    let _name = args[0];
    let _group = args[1];
    let style = FixStyle::from(args[2]);

    let rest = &args[3..];

    match style {
        FixStyle::Nve => parse_no_args(rest),
        FixStyle::Nvt => parse_nh_fixes(rest, style),

        _ => todo!("Parsing for this fix style is not yet implemented."),
    }
}

/// Generic Parsing of the Nose-Hoover Fixes
/// TODO Finish ME
fn parse_nh_fixes(_rest: &[&str], _fix_style: FixStyle) -> Result<(), InvalidArguments> {
    todo!()
}

/// Parse a fix that does not expect any arguments
fn parse_no_args(rest: &[&str]) -> Result<(), InvalidArguments> {
    if !rest.is_empty() {
        Err(InvalidArguments::IncorrectNumberArguments {
            provided: rest.len() + 3,
            expected: 3,
        })
    } else {
        Ok(())
    }
}

/// Acceptable arguments for LAMMPS commands
/// TODO Convert from a tree-sitter nodes to avoid reparsing over again!!!
pub enum Argument {
    /// Expression from LAMMPS
    Expression,
    String,
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    /// TODO make this hold and `Ident` struct, from `crate::identifinder` module
    UnderscoreIdent,
    /// Variables within curly braces
    /// TODO make this hold and `Ident` struct, from `crate::identifinder` module
    CurlyVar,
}

// #[derive(Parser)]
// struct FixNVT {
//     // Required Args
//     fix_id: String,
//     group_id: String,
//     #[arg(long)]
//     temp: Vec<f64>,
//     /// Start stop damp
//     // Optional Args
//     tchain: u64,
//     mtk: bool,
//     drag: f64,
//     update: DipoleUpdate,
// }

// #[derive(Clone)]
// enum DipoleUpdate {
//     Dipole,
//     DipoleDLM,
// }

// impl From<&str> for DipoleUpdate {
//     fn from(s: &str) -> Self {
//         match s {
//             "dipole" => Self::Dipole,
//             "dipole/dlm" => Self::DipoleDLM,
//             _ => panic!("Unknown dipole update type"),
//         }
//     }
// }

#[cfg(test)]
mod tests {
    use tree_sitter::{Query, QueryCursor};

    #[test]
    fn read_ts_node() {
        use tree_sitter::Parser;
        let mut parser = Parser::new();
        parser.set_language(tree_sitter_lammps::language()).unwrap();
        let source = include_bytes!("../just_a_fix.lmp");
        let tree = parser.parse(source, None).unwrap();
        let mut tree_cursor = tree.walk();
        // Getting the arguments node

        let query = Query::new(tree_sitter_lammps::language(), "(args) @args").unwrap();

        let mut cursor = QueryCursor::new();
        let matches = cursor.captures(&query, tree.root_node(), source.as_slice());

        for mtch in matches {
            dbg!(&mtch.0.captures[0]
                .node
                .children(&mut tree_cursor)
                .collect::<Vec<_>>());
        }

        panic!("Failing test on purpose!!!")
    }
}

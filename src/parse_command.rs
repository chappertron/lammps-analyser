// Parse LAMMPS commands into a struct using Clap?
// TODO I am not sure this will work. Primarily because there
// are no -- flags. I think I will need to use a custom parser :(

// use clap::Parser;

use thiserror::Error;

use crate::ast::{Argument, FixDef};
use crate::fix_styles::FixStyle;

#[derive(Debug, Error, PartialEq, Eq)]
pub enum InvalidArguments {
    #[error("Incorrect number of arguments: {provided} expected {expected}")]
    IncorrectNumberArguments { provided: usize, expected: usize },
    #[error(
        "Incorrect keyword arguments for `{kwarg}`. 
         Expected {n_expected} arguments: {expected}. 
         Only {n_provided} provided."
    )]
    MissingKwargField {
        kwarg: String,
        expected: String,
        n_expected: usize,
        n_provided: usize,
    },
    #[error("Incorrect argument type. Expected: {expected}. Provided: {provided}.")]
    IncorrectType { expected: String, provided: String },
    #[error("{0}")]
    /// Error for specific fixes
    Custom(String),
}

/// Parse
pub fn parse_fix(fix: FixDef) -> Result<(), InvalidArguments> {
    // if args.len() < 3 {
    //     return Err(InvalidArguments::IncorrectNumberArguments {
    //         provided: args.len(),
    //         expected: 3,
    //     });
    // }

    // TODO use these
    let style = fix.fix_style;

    match style {
        FixStyle::Nve => parse_no_args(&fix),
        FixStyle::Nvt => parse_nh_fixes(&fix),

        _ => todo!("Parsing for this fix style is not yet implemented."),
    }
}

/// Generic Parsing of the Nose-Hoover Fixes
/// TODO Finish ME
fn parse_nh_fixes(fix: &FixDef) -> Result<(), InvalidArguments> {
    let args = &fix.args;

    // Iterate through and check validity of the argument
    let mut iter = args.iter();

    // Try either the LAMMPS way or to use peek with a while let loop?

    while let Some(arg) = iter.next() {
        match arg {
            Argument::ArgName(kwarg) => {
                // TODO Convert into match block
                if kwarg == "temp" {
                    // TODO check if there are 3 more elements
                    kwarg_expected_floats(&mut iter, kwarg, 3, "<Tstart> <Tstop> <Tdamp>")?;
                } else if matches!(
                    kwarg.as_ref(),
                    "iso" | "aniso" | "tri" | "x" | "y" | "z" | "xy" | "xz" | "yz"
                ) {
                    if !matches!(fix.fix_style, FixStyle::Npt | FixStyle::Nph) {
                        Err(InvalidArguments::Custom(format!(
                            "{} is not a valid keyword for fix {}",
                            kwarg, fix.fix_style
                        )))?
                    }
                    kwarg_expected_floats(&mut iter, kwarg, 3, "<Pstart> <Pstop> <Pdamp>")?;
                } else {
                    todo!("Other keyword arguments are not yet implemented.")
                }
            }

            _ => todo!(),
        }
    }

    Ok(())
}

fn kwarg_expected_floats<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), InvalidArguments> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x {
                Argument::Float(_) => (),
                Argument::Int(_) => (),
                Argument::VarRound(_) => (),
                Argument::VarCurly(_) => (),
                _ => Err(InvalidArguments::IncorrectType {
                    expected: "float".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(InvalidArguments::MissingKwargField {
                kwarg: kwarg.into(),
                n_expected,
                n_provided: i,
                expected: expected_args.into(),
            })?
        }
    }

    Ok(())
}

// fn kwarg_expected_float_expr(
//     &mut iter: impl Iterator,
//     kwarg: &str,
//     n_expected: usize,
// ) -> Result<(), InvalidArguments> {
// }

/// Parse a fix that does not expect any arguments
pub fn parse_no_args(fix: &FixDef) -> Result<(), InvalidArguments> {
    if !fix.args.is_empty() {
        Err(InvalidArguments::IncorrectNumberArguments {
            provided: fix.args.len() + 3,
            expected: 3,
        })
    } else {
        Ok(())
    }
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
    use super::*;
    use crate::fix_styles::FixStyle;
    use tree_sitter::{Parser, Query, QueryCursor};

    fn setup_parser() -> Parser {
        let mut parser = Parser::new();
        parser.set_language(tree_sitter_lammps::language()).unwrap();
        parser
    }
    #[test]
    fn valid_zero_arg_fix() {
        let mut parser = setup_parser();
        let text = "fix NVE all nve";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes());

        assert_eq!(fix.fix_id.name, "NVE");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Nve);
        assert!(fix.args.is_empty());

        assert_eq!(parse_no_args(&fix), Ok(()));
    }

    #[test]
    fn invalid_zero_arg_fix() {
        let mut parser = setup_parser();
        let text = "fix NVE all nve asdfas";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes());

        assert_eq!(fix.fix_id.name, "NVE");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Nve);
        assert!(!fix.args.is_empty());

        assert_eq!(
            parse_no_args(&fix),
            Err(InvalidArguments::IncorrectNumberArguments {
                provided: 4,
                expected: 3
            })
        );
    }

    #[test]
    fn valid_nvt() {
        let mut parser = setup_parser();
        let text = "fix NVT all nvt temp 1 $(v_T*1.5) $(100*dt)";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes());

        assert_eq!(fix.fix_id.name, "NVT");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Nvt);
        assert!(!fix.args.is_empty());

        assert_eq!(parse_nh_fixes(&fix), Ok(()));
    }

    #[test]
    fn valid_npt() {
        let mut parser = setup_parser();
        let text = "fix NPT all npt temp 1 $(v_T*1.5) $(100*dt) iso 1 $(v_p*1.5) ${pdamp}";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes());

        assert_eq!(fix.fix_id.name, "NPT");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Npt);
        assert!(!fix.args.is_empty());

        assert_eq!(dbg!(parse_nh_fixes(&fix)), Ok(()));
    }

    #[test]
    fn invalid_nvt() {
        let mut parser = setup_parser();
        let text = "fix NVT all nvt temp 1.0 $(v_T*1.5) TEMP";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes());

        assert_eq!(fix.fix_id.name, "NVT");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Nvt);
        assert!(!fix.args.is_empty());

        assert_eq!(parse_nh_fixes(&fix), Ok(()));
    }

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

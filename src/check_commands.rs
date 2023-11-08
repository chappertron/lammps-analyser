// Checks the validity of fix arguments

use lsp_types::DiagnosticSeverity;
use owo_colors::OwoColorize;
use thiserror::Error;
use tree_sitter::{Point, Range};

use crate::ast::{Argument, FixDef};
use crate::diagnostic_report::ReportSimple;
use crate::fix_styles::FixStyle;
use crate::utils::point_to_position;

#[derive(Debug, Error, PartialEq, Eq)]
#[error("{}:{}: {err_type} for fix: {fix_style}",range.start_point.row+1,range.start_point.column+1)]
pub struct InvalidArguments {
    pub(crate) err_type: InvalidArgumentsType,
    pub(crate) range: Range,
    pub(crate) fix_style: FixStyle,
}

impl InvalidArguments {
    pub fn new(err_type: InvalidArgumentsType, range: Range, fix_style: FixStyle) -> Self {
        InvalidArguments {
            err_type,
            range,
            fix_style,
        }
    }

    pub fn range(&self) -> Range {
        self.range
    }
    pub fn start(&self) -> Point {
        self.range.start_point
    }
    pub fn end(&self) -> Point {
        self.range.end_point
    }
}

impl From<InvalidArguments> for lsp_types::Diagnostic {
    fn from(value: InvalidArguments) -> Self {
        lsp_types::Diagnostic {
            range: lsp_types::Range {
                start: point_to_position(&(value.start())),
                end: point_to_position(&(value.end())),
            },
            severity: Some(DiagnosticSeverity::ERROR),
            message: value.to_string(),

            ..Default::default()
        }
    }
}

impl ReportSimple for InvalidArguments {
    fn make_simple_report(&self) -> String {
        format!(
            "{}:{}: {} for fix: {}",
            self.start().row + 1,
            self.start().column + 1,
            self.err_type.bright_red(),
            self.fix_style.bright_red(),
        )
    }
}

#[derive(Debug, Error, PartialEq, Eq)]
pub enum InvalidArgumentsType {
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

    #[error("{} is not a valid keyword for fix {}", kwarg, fix_style)]
    InvalidKeyword { kwarg: String, fix_style: FixStyle },
    #[error("Invalid option `{}` for keyword `{}`. Valid options: [{}].", kwarg, provided,options.join(","))]
    InvalidOption {
        kwarg: String,
        provided: String,
        options: Vec<String>,
    },
    #[error("Invalid style for fix: {0}")]
    InvalidStyle(String),
}

/// Check the arguments of a fix
pub fn check_fix(fix: &FixDef) -> Result<(), InvalidArguments> {
    let style = fix.fix_style;

    // TODO: Check group of the fix

    match style {
        //  TODO: Report the invalid styles name

        // FixStyle::InvalidFixStyle => Err(InvalidArguments {
        //     err_type: InvalidArgumentsType::InvalidStyle(fix.fix_style.to_string()),
        //     range: fix.range(),
        //     fix_style: fix.fix_style,
        // }),

        // `Ok`, so duplicates aren't raised
        FixStyle::InvalidFixStyle => Ok(()),

        FixStyle::Nve => parse_no_args(fix).map_err(|x| InvalidArguments {
            err_type: x,
            range: fix.range(),
            fix_style: style,
        }),
        // TODO: See if other fix styles share the same arguments
        FixStyle::Nvt | FixStyle::Npt | FixStyle::Nph => {
            parse_nh_fixes(fix).map_err(|x| InvalidArguments {
                err_type: x,
                range: fix.range(),
                fix_style: style,
            })
        }
        // Fallback to checking only positional arguments
        style => check_n_positional(fix, style.n_positional_args()).map_err(|x| InvalidArguments {
            err_type: x,
            range: fix.range(),
            fix_style: style,
        }),
        // _ => Ok(()), // Ignore unchecked fixes
        // _ => Err(InvalidArguments::Custom(
        //     "Parsing for this fix style is not yet implemented.".into(),
        // )),
    }
}

/// Parse Fix with at least n positional arguments
fn check_n_positional(fix: &FixDef, n_args: usize) -> Result<(), InvalidArgumentsType> {
    if fix.args.len() < n_args {
        Err(InvalidArgumentsType::IncorrectNumberArguments {
            provided: fix.args.len(),
            expected: n_args,
        })
    } else {
        Ok(())
    }
}

/// Generic Parsing of the Nose-Hoover Fixes
/// TODO: Finish ME
fn parse_nh_fixes(fix: &FixDef) -> Result<(), InvalidArgumentsType> {
    let args = &fix.args;

    // Iterate through and check validity of the argument
    let mut iter = args.iter();

    let is_barostatting = matches!(fix.fix_style, FixStyle::Npt | FixStyle::Nph);
    // Try either the LAMMPS way or to use peek with a while let loop?
    let barostat_only = |kwarg: &str| {
        if !is_barostatting {
            Err(InvalidArgumentsType::InvalidKeyword {
                kwarg: kwarg.to_string(),
                fix_style: fix.fix_style,
            })
        } else {
            Ok(())
        }
    };
    while let Some(arg) = iter.next() {
        match arg {
            Argument::ArgName(kwarg) if kwarg == "temp" => {
                // TODO: check if there are 3 more elements
                kwarg_expected_floats(&mut iter, kwarg, 3, "<Tstart> <Tstop> <Tdamp>")?;
            }
            Argument::ArgName(kwarg)
                if matches!(
                    kwarg.as_ref(),
                    "iso" | "aniso" | "tri" | "x" | "y" | "z" | "xy" | "xz" | "yz"
                ) =>
            {
                kwarg_expected_floats(&mut iter, kwarg, 3, "<Pstart> <Pstop> <Pdamp>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "couple" => {
                barostat_only(kwarg)?;
                kwarg_expected_enum(&mut iter, kwarg, 1, &["none", "xyz", "xy", "xz", "yz"])?;
            }
            Argument::ArgName(kwarg) if kwarg == "tchain" => {
                kwarg_expected_floats(&mut iter, kwarg, 1, "<N> chain.")?;
            }
            Argument::ArgName(kwarg) if kwarg == "pchain" => {
                barostat_only(kwarg)?;
                kwarg_expected_floats(&mut iter, kwarg, 1, "<N> chain.")?;
            }
            Argument::ArgName(kwarg) if kwarg == "mtk" => {
                barostat_only(kwarg)?;
                kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "tloop" => {
                kwarg_expected_floats(&mut iter, kwarg, 1, "<N> sub-cycles")?;
            }
            Argument::ArgName(kwarg) if kwarg == "ploop" => {
                barostat_only(kwarg)?;
                kwarg_expected_floats(&mut iter, kwarg, 1, "<N> sub-cycles")?;
            }
            Argument::ArgName(kwarg) if kwarg == "nreset" => {
                barostat_only(kwarg)?;
                kwarg_expected_floats(&mut iter, kwarg, 1, "<N> reset")?;
            }
            Argument::ArgName(kwarg) if kwarg == "drag" => {
                barostat_only(kwarg)?;
                kwarg_expected_floats(&mut iter, kwarg, 1, "<Df> drag factor")?;
            }
            Argument::ArgName(kwarg) if kwarg == "ptemp" => {
                barostat_only(kwarg)?;
                kwarg_expected_floats(&mut iter, kwarg, 1, "<Ttarget>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "dilate" => {
                barostat_only(kwarg)?;
                // TODO: convert to a group
                kwarg_expected_str(&mut iter, kwarg, 1, "<dilate-group-ID>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "scalexy" => {
                barostat_only(kwarg)?;
                kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "scaleyz" => {
                barostat_only(kwarg)?;
                kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "scalexz" => {
                barostat_only(kwarg)?;
                kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }

            Argument::ArgName(kwarg) if kwarg == "flip" => {
                barostat_only(kwarg)?;
                kwarg_expected_bool(&mut iter, kwarg, 1, "<value>")?;
            }

            Argument::ArgName(kwarg) if kwarg == "fixedpoint" => {
                barostat_only(kwarg)?;
                kwarg_expected_floats(&mut iter, kwarg, 1, "<x> <y> <z>")?;
            }
            Argument::ArgName(kwarg) if kwarg == "update" => {
                kwarg_expected_enum(&mut iter, kwarg, 1, &["dipole", "dipole/dlm"])?
            }

            Argument::ArgName(kwarg) => Err(InvalidArgumentsType::Custom(format!(
                "Unknown kwarg argument: {kwarg}",
            )))?,

            _ => Err(InvalidArgumentsType::Custom(format!(
                "Unknown argument: {arg}",
            )))?,
        }
    }

    Ok(())
}

/// Parse n keywords that match in the provided slice
fn kwarg_expected_enum<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_args: usize,
    options: &[&str],
) -> Result<(), InvalidArgumentsType> {
    for i in 0..n_args {
        if let Some(x) = iter.next() {
            match x {
                Argument::ArgName(x) => {
                    if !options.contains(&x.as_ref()) {
                        Err(InvalidArgumentsType::InvalidOption {
                            kwarg: kwarg.into(),
                            provided: x.to_string(),
                            options: options.iter().map(|&x| x.to_string()).collect(),
                        })?
                    }
                }
                _ => Err(InvalidArgumentsType::IncorrectType {
                    expected: "string-like".into(),
                    provided: x.to_string(),
                })?,
            }
        } else {
            Err(InvalidArgumentsType::MissingKwargField {
                kwarg: kwarg.into(),
                expected: format!("One of: {}.", options.join(", ")),
                n_expected: n_args,
                n_provided: i,
            })?
        }
    }
    Ok(())
}

fn kwarg_expected_floats<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), InvalidArgumentsType> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x {
                Argument::Float(_) => (),
                Argument::Int(_) => (),
                Argument::VarRound(_) => (),
                Argument::VarCurly(_) => (),
                _ => Err(InvalidArgumentsType::IncorrectType {
                    expected: "float".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(InvalidArgumentsType::MissingKwargField {
                kwarg: kwarg.into(),
                n_expected,
                n_provided: i,
                expected: expected_args.into(),
            })?
        }
    }

    Ok(())
}

fn kwarg_expected_bool<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), InvalidArgumentsType> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x {
                Argument::Bool(_) => (),
                _ => Err(InvalidArgumentsType::IncorrectType {
                    expected: "bool".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(InvalidArgumentsType::MissingKwargField {
                kwarg: kwarg.into(),
                n_expected,
                n_provided: i,
                expected: expected_args.into(),
            })?
        }
    }

    Ok(())
}

fn kwarg_expected_str<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), InvalidArgumentsType> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x {
                Argument::ArgName(_) => (),
                _ => Err(InvalidArgumentsType::IncorrectType {
                    expected: "string".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(InvalidArgumentsType::MissingKwargField {
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
pub fn parse_no_args(fix: &FixDef) -> Result<(), InvalidArgumentsType> {
    if !fix.args.is_empty() {
        Err(InvalidArgumentsType::IncorrectNumberArguments {
            provided: fix.args.len(),
            expected: 0,
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
    use crate::ast::from_node::FromNode;
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
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

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
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "NVE");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Nve);
        assert!(!fix.args.is_empty());

        assert_eq!(
            parse_no_args(&fix),
            Err(InvalidArgumentsType::IncorrectNumberArguments {
                provided: 1,
                expected: 0
            })
        );
    }

    #[test]
    fn valid_nvt() {
        let mut parser = setup_parser();
        let text = "fix NVT all nvt temp 1 $(v_T*1.5) $(100*dt)";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

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
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "NPT");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Npt);
        assert!(!fix.args.is_empty());

        assert_eq!(dbg!(parse_nh_fixes(&fix)), Ok(()));
    }

    #[test]
    fn complicated_nph() {
        let mut parser = setup_parser();
        let text = "fix 2 ice nph x 1.0 1.0 0.5 y 2.0 2.0 0.5 z 3.0 3.0 0.5 yz 0.1 0.1 0.5 xz 0.2 0.2 0.5 xy 0.3 0.3 0.5 nreset 1000";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "2");
        assert_eq!(fix.group_id, "ice");
        assert_eq!(fix.fix_style, FixStyle::Nph);
        assert!(!fix.args.is_empty());

        assert_eq!(dbg!(parse_nh_fixes(&fix)), Ok(()));
    }

    #[test]
    fn invalid_nvt() {
        let mut parser = setup_parser();
        let text = "fix NVT all nvt temp 1.0 $(v_T*1.5) TEMP";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert_eq!(fix.fix_id.name, "NVT");
        assert_eq!(fix.group_id, "all");
        assert_eq!(fix.fix_style, FixStyle::Nvt);
        assert!(!fix.args.is_empty());

        // TODO: Turn into an error...
        assert!(parse_nh_fixes(&fix).is_err());
    }

    #[test]
    fn check_n_args() {
        let mut parser = setup_parser();
        let text = "fix profw water  ave/chunk 1 $(v_nsteps) $(v_nsteps) lwater  density/number density/mass  temp file ${outputname}water.profiles adof 2";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        dbg!(&fix);

        assert!(check_fix(&fix).is_ok())
    }

    #[test]
    fn bad_check_n_args() {
        let mut parser = setup_parser();
        let text = "fix profw water  ave/chunk 1 $(v_nsteps) $(v_nsteps)";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert!(dbg!(check_fix(&fix)).is_err())
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

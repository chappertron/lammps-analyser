//! Utilities for checking the arguments of commands/fixes/computes.
//!
//!
//!
//!
//!

use crate::ast::{Argument, ArgumentKind, FixDef};

use super::invalid_arguments;

/// Check that the provided arguments are one of several provided options.
///
/// Takes the `n_args` arguments from the iterator and checks each one is in `options`
pub(crate) fn kwarg_expected_enum<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &'static str,
    n_args: usize,
    options: &[&'static str],
) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    for i in 0..n_args {
        if let Some(x) = iter.next() {
            match &x.kind {
                ArgumentKind::ArgName(x) | ArgumentKind::Word(x) => {
                    if !options.contains(&x.as_ref()) {
                        Err(invalid_arguments::InvalidArgumentsType::InvalidOption {
                            kwarg,
                            provided: x.to_string(),
                            options: options.to_vec(),
                        })?
                    }
                }
                _ => Err(invalid_arguments::InvalidArgumentsType::IncorrectType {
                    expected: "string-like".into(),
                    provided: x.to_string(),
                })?,
            }
        } else {
            Err(invalid_arguments::InvalidArgumentsType::MissingKwargField {
                kwarg: kwarg.into(),
                expected: format!("One of: {}.", options.join(", ")),
                n_expected: n_args,
                n_provided: i,
            })?
        }
    }
    Ok(())
}

pub(crate) fn kwarg_expected_floats<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x.kind {
                ArgumentKind::Float(_) => (),
                ArgumentKind::Int(_) => (),
                ArgumentKind::VarRound(_) => (),
                ArgumentKind::VarCurly(_) => (),
                ArgumentKind::SimpleExpansion(_) => (),
                _ => Err(invalid_arguments::InvalidArgumentsType::IncorrectType {
                    expected: "float".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(invalid_arguments::InvalidArgumentsType::MissingKwargField {
                kwarg: kwarg.into(),
                n_expected,
                n_provided: i,
                expected: expected_args.into(),
            })?
        }
    }

    Ok(())
}

pub(crate) fn kwarg_expected_bool<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x.kind {
                ArgumentKind::Bool(_) => (),
                _ => Err(invalid_arguments::InvalidArgumentsType::IncorrectType {
                    expected: "bool".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(invalid_arguments::InvalidArgumentsType::MissingKwargField {
                kwarg: kwarg.into(),
                n_expected,
                n_provided: i,
                expected: expected_args.into(),
            })?
        }
    }

    Ok(())
}

pub(crate) fn kwarg_expected_str<'a>(
    iter: &mut impl Iterator<Item = &'a Argument>,
    kwarg: &str,
    n_expected: usize,
    expected_args: &str,
) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    for i in 0..n_expected {
        if let Some(x) = iter.next() {
            match x.kind {
                ArgumentKind::ArgName(_) | ArgumentKind::Word(_) => (),
                _ => Err(invalid_arguments::InvalidArgumentsType::IncorrectType {
                    expected: "string".into(),
                    provided: x.to_string(),
                })?,
            };
        } else {
            Err(invalid_arguments::InvalidArgumentsType::MissingKwargField {
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
pub fn parse_no_args(fix: &FixDef) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    if !fix.args.is_empty() {
        Err(
            invalid_arguments::InvalidArgumentsType::IncorrectNumberArguments {
                provided: fix.args.len(),
                expected: 0,
            },
        )
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ast::from_node::FromNode, check_commands::utils::parse_no_args, styles::FixStyle,
        utils::testing::parse,
    };

    #[test]
    fn valid_zero_arg_fix() {
        let text = "fix NVE all nve";
        let tree = parse(text);
        let node = tree.root_node().child(0).unwrap();
        let fix = FixDef::from_node(&node, text).unwrap();

        assert_eq!(fix.fix_id.name, "NVE");
        assert_eq!(fix.group_id.contents, "all");
        assert_eq!(fix.fix_style, FixStyle::Nve);
        assert!(fix.args.is_empty());

        assert_eq!(parse_no_args(&fix), Ok(()));
    }

    #[test]
    fn invalid_zero_arg_fix() {
        let text = "fix NVE all nve asdfas";
        let tree = parse(text);
        let node = tree.root_node().child(0).unwrap();
        let fix = FixDef::from_node(&node, text).unwrap();

        assert_eq!(fix.fix_id.name, "NVE");
        assert_eq!(fix.group_id.contents, "all");
        assert_eq!(fix.fix_style, FixStyle::Nve);
        assert!(!fix.args.is_empty());

        assert_eq!(
            parse_no_args(&fix),
            Err(
                invalid_arguments::InvalidArgumentsType::IncorrectNumberArguments {
                    provided: 1,
                    expected: 0
                }
            )
        );
    }
}

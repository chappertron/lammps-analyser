use crate::ast::FixDef;
use crate::commands::CommandName;
use crate::compute_styles::ComputeStyle;

use super::invalid_arguments;
use super::utils::parse_no_args;

use crate::fix_styles::FixStyle;

/// enumerate between styles and commands.

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CommandAndStyle {
    /// A compute style
    ComputeStyle(ComputeStyle),
    /// A fix style
    FixStyle(FixStyle),
    /// A commandName
    CommandName(CommandName),
}

impl std::fmt::Display for CommandAndStyle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ComputeStyle(compute) => write![f, "compute `{}`", compute],
            Self::FixStyle(fix) => write![f, "fix `{}`", fix],
            Self::CommandName(command) => write![f, "command `{}`", command],
        }
    }
}

/// Check the arguments of a fix
pub fn check_fix(fix: &FixDef) -> Result<(), invalid_arguments::InvalidArguments> {
    let style = fix.fix_style;

    // TODO: Check group of the fix. This would likely be a separate function or require more
    // arguments in this function.

    match style {
        //  TODO: Report the invalid styles name here
        // `Ok`, so duplicates aren't raised. Also checked elsewhere...
        FixStyle::InvalidStyle => Ok(()),

        FixStyle::Nve => parse_no_args(fix).map_err(|x| invalid_arguments::InvalidArguments {
            err_type: x,
            range: fix.range(),
            style: CommandAndStyle::FixStyle(style),
        }),
        // TODO: See if other fix styles share the same arguments
        FixStyle::Nvt | FixStyle::Npt | FixStyle::Nph => {
            nose_hoover::parse_nh_fixes(fix).map_err(|x| invalid_arguments::InvalidArguments {
                err_type: x,
                range: fix.range(),
                style: CommandAndStyle::FixStyle(style),
            })
        }
        // Fallback to checking only positional arguments
        style => check_n_positional(fix, style.n_positional_args()).map_err(|x| {
            invalid_arguments::InvalidArguments {
                err_type: x,
                range: fix.range(),
                style: CommandAndStyle::FixStyle(style),
            }
        }),
    }
}

// NOTE: Tests are in the appropiate sub modules for the given fix.

/// Parse Fix with at least n positional arguments
pub(crate) fn check_n_positional(
    fix: &FixDef,
    n_args: usize,
) -> Result<(), invalid_arguments::InvalidArgumentsType> {
    if fix.args.len() < n_args {
        Err(
            invalid_arguments::InvalidArgumentsType::IncorrectNumberArguments {
                provided: fix.args.len(),
                expected: n_args,
            },
        )
    } else {
        Ok(())
    }
}

pub(crate) mod nose_hoover;

#[cfg(test)]
mod tests {
    use crate::{ast::from_node::FromNode, check_commands::fixes, utils::testing::setup_parser};

    use super::*;

    #[test]
    fn check_n_args() {
        let mut parser = setup_parser();
        let text = "fix profw water  ave/chunk 1 $(v_nsteps) $(v_nsteps) lwater  density/number density/mass  temp file ${outputname}water.profiles adof 2";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        dbg!(&fix);

        assert!(fixes::check_fix(&fix).is_ok())
    }

    #[test]
    fn bad_check_n_args() {
        let mut parser = setup_parser();
        let text = "fix profw water  ave/chunk 1 $(v_nsteps) $(v_nsteps)";
        let tree = parser.parse(text, None).unwrap();
        let node = tree.root_node().child(0).unwrap().child(0).unwrap();
        let fix = FixDef::from_node(&node, text.as_bytes()).unwrap();

        assert!(dbg!(fixes::check_fix(&fix)).is_err())
    }
}

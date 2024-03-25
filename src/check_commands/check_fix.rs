use super::check_n_positional;
use super::invalid_arguments;

use super::parse_nh_fixes;

use super::parse_no_args;

use crate::fix_styles::FixStyle;

use crate::ast::FixDef;

/// Check the arguments of a fix
pub fn check_fix(fix: &FixDef) -> Result<(), invalid_arguments::InvalidArguments> {
    let style = fix.fix_style;

    // TODO: Check group of the fix. This qould likely be a separate function

    match style {
        //  TODO: Report the invalid styles name

        // FixStyle::InvalidFixStyle => Err(InvalidArguments {
        //     err_type: InvalidArgumentsType::InvalidStyle(fix.fix_style.to_string()),
        //     range: fix.range(),
        //     fix_style: fix.fix_style,
        // }),

        // `Ok`, so duplicates aren't raised
        FixStyle::InvalidFixStyle => Ok(()),

        FixStyle::Nve => parse_no_args(fix).map_err(|x| invalid_arguments::InvalidArguments {
            err_type: x,
            range: fix.range(),
            fix_style: style,
        }),
        // TODO: See if other fix styles share the same arguments
        FixStyle::Nvt | FixStyle::Npt | FixStyle::Nph => {
            parse_nh_fixes(fix).map_err(|x| invalid_arguments::InvalidArguments {
                err_type: x,
                range: fix.range(),
                fix_style: style,
            })
        }
        // Fallback to checking only positional arguments
        style => check_n_positional(fix, style.n_positional_args()).map_err(|x| {
            invalid_arguments::InvalidArguments {
                err_type: x,
                range: fix.range(),
                fix_style: style,
            }
        }),
        // _ => Ok(()), // Ignore unchecked fixes
        // _ => Err(InvalidArguments::Custom(
        //     "Parsing for this fix style is not yet implemented.".into(),
        // )),
    }
}

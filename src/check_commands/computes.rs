use crate::ast::ComputeDef;

use super::{fixes::CommandAndStyle, invalid_arguments::InvalidArguments};

/// Checks whether the provided compute definition has valid arguments.
///
/// All computes are checked for a minimum number of arguments they are expected to have, as
/// defined in [`crate::compute_styles::ComputeStyle`]
///
/// Returns `Ok(())` if there are no issues, otherwise returns the appropiate error.
pub fn check_compute(compute: &ComputeDef) -> Result<(), InvalidArguments> {
    if compute.args.len() < compute.compute_style.n_positional_args() {
        return Err(InvalidArguments::new(
            super::invalid_arguments::InvalidArgumentsType::IncorrectNumberArguments {
                provided: compute.args.len(),
                expected: compute.compute_style.n_positional_args(),
            },
            compute.range(),
            CommandAndStyle::ComputeStyle(compute.compute_style),
        ));
    }

    Ok(())
}

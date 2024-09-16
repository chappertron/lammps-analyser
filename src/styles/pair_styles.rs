#![allow(unused)]

include!(concat!(env!("OUT_DIR"), "/styles/pair_styles.rs"));

/// Containing information about the PairStyle
pub struct PairStyleInfo {
    /// The minimum number of arguments the `pair_style` command takes
    pub min_args: usize,
    /// The maximum Number of arguments the `pair_style` command takes
    pub max_args: Option<usize>,
    /// The minimum number of coefficients that need to be provided to the `pair_coeff`
    /// command.
    pub min_coeffs: usize,
    /// The maximum number of coefficients that need to be provided to the `pair_coeff`
    /// command.
    pub max_coeffs: Option<usize>,
    pub arg_context: &'static str,
    pub coeff_context: &'static str,
}

impl PairStyle {
    // TODO: Write a JSON file or somehow generate this?

    /// Information about the style
    pub fn info(&self) -> Option<PairStyleInfo> {
        use PairStyle as PS;
        match self {
            PS::LjCut => Some(PairStyleInfo {
                min_args: 1,
                max_args: Some(1),
                min_coeffs: 2,
                max_coeffs: Some(3),
                arg_context: "cutoff",
                coeff_context: "epsilon sigma (cutoff)",
            }),
            PS::LjCutCoulCut => Some(PairStyleInfo {
                min_args: 1,
                max_args: Some(2),
                min_coeffs: 2,
                // Optional extra cutoffs
                max_coeffs: Some(4),
                arg_context: "cutoff (cutoff2)",
                coeff_context: "epsilon sigma (cutoff) (cutoff2)",
            }),
            PS::LjCutCoulLong | PS::LjCutCoulMsm => Some(PairStyleInfo {
                min_args: 1,
                max_args: Some(2),
                min_coeffs: 2,
                // Optional extra cutoffs. Can't specify coul cutoff seperately for these
                // styles
                max_coeffs: Some(3),
                arg_context: "cutoff (cutoff2)",
                coeff_context: "epsilon sigma (cutoff)",
            }),
            PS::LjCutCoulDebye | PS::LjCutCoulDsf | PS::LjCutCoulWolf => Some(PairStyleInfo {
                // alpha/kappa cut1 (cut2)
                min_args: 2,
                max_args: Some(3),
                min_coeffs: 2,
                max_coeffs: Some(4),
                arg_context: "alpha/kappa cutoff (cutoff2)",
                coeff_context: "epsilon sigma (cutoff) (cutoff2)",
            }),
            PS::LjCutTip4pCut | PS::LjCutTip4pLong => Some(PairStyleInfo {
                // otype htype btype atype qdist cutoff (cutoff2)
                min_args: 6,
                max_args: Some(7),
                // sigma, epsilon (cut_lj)
                min_coeffs: 2,
                // Optional pairwise lj cut-off
                max_coeffs: Some(3),
                arg_context: "otype htype btype atype qdist cutoff (cutoff2)",
                coeff_context: "epsilon sigma (lj_cutoff)",
            }),

            PS::Hybrid => Some(PairStyleInfo {
                min_args: 1,
                // Can in principle have any number of args
                max_args: None,
                min_coeffs: 1,
                // Can in principle have any number of coeffs
                max_coeffs: None,
                arg_context: "style1 args style2 args...",
                coeff_context: "style args",
            }),

            _ => None,
        }
    }

    pub fn min_args(&self) -> Option<usize> {
        self.info().map(|i| i.min_args)
    }

    pub fn max_args(&self) -> Option<usize> {
        self.info().map(|i| i.max_args)?
    }

    pub fn min_coeffs(&self) -> Option<usize> {
        self.info().map(|i| i.min_coeffs)
    }

    pub fn max_coeffs(&self) -> Option<usize> {
        self.info().map(|i| i.max_coeffs)?
    }

    pub fn arg_context(&self) -> Option<&str> {
        self.info().map(|i| i.arg_context)
    }

    pub fn coeff_context(&self) -> Option<&str> {
        self.info().map(|i| i.coeff_context)
    }
}

#[cfg(test)]
mod tests {

    #[test]
    #[ignore = "test is incomplete."]
    fn correct_number_of_args_in_context() {
        todo!("try all styles that have infos and check if the context message matches the expected number of args and funcs.

");
    }
}

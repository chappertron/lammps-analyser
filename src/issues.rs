/// Types for LAMMPS issues and errorsmm
use crate::lammps_errors::{LammpsError, Warnings};
use owo_colors::OwoColorize;
use thiserror::Error;

#[derive(Error, Debug, Clone, PartialEq, Eq)]
pub enum Issue {
    /// Error in LAMMPS input file
    #[error("Error: {0}")]
    Error(LammpsError),
    #[error("Warning:")]
    Warning(Warnings),
    /// TODO: NOT CURRENTLY IMPLEMENTED
    #[error("Info:")]
    Info,
}

impl From<Warnings> for Issue {
    fn from(v: Warnings) -> Self {
        Self::Warning(v)
    }
}

impl From<LammpsError> for Issue {
    fn from(v: LammpsError) -> Self {
        Self::Error(v)
    }
}

impl From<Issue> for lsp_types::Diagnostic {
    fn from(value: Issue) -> Self {
        match value {
            Issue::Error(e) => e.into(),
            _ => unimplemented!("Warning and Info not implemented yet"),
        }
    }
}

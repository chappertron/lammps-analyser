/// Types for LAMMPS issues and errorsmm
use crate::{diagnostic_report::ReportSimple, lammps_errors::LammpsError};
use owo_colors::OwoColorize;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum Issue {
    /// Error in LAMMPS input file
    #[error("Error: {0}")]
    Error(LammpsError),
    /// TODO NOT CURRENTLY IMPLEMENTED
    #[error("Warning:")]
    Warning,
    /// TODO: NOT CURRENTLY IMPLEMENTED
    #[error("Info:")]
    Info,
}

impl From<LammpsError> for Issue {
    fn from(v: LammpsError) -> Self {
        Self::Error(v)
    }
}

impl ReportSimple for Issue {
    fn make_simple_report(&self) -> String {
        match self {
            Issue::Error(e) => format!("{} {}", "Error:".bright_red(), e.make_simple_report()),
            Issue::Warning => "Warning".to_string(),
            Issue::Info => "Info".to_string(),
        }
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

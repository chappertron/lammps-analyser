use crate::check_commands::InvalidArguments;
use crate::check_styles::InvalidStyle;
use crate::diagnostic_report::ReportSimple;
use crate::error_finder::SyntaxError;
use crate::identifinder::{UndefinedIdent, UnusedIdent};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum LammpsError {
    #[error("{0}")]
    SyntaxError(SyntaxError),
    #[error("{0}")]
    InvalidStyle(InvalidStyle),
    #[error("{0}")]
    UndefinedIdent(UndefinedIdent),
    #[error("{0}")]
    InvalidArguments(InvalidArguments),
}

#[derive(Error, Debug)]
pub enum Warnings {
    #[error("{0}")]
    UnusedIdent(UnusedIdent),
}

impl From<UnusedIdent> for Warnings {
    fn from(v: UnusedIdent) -> Self {
        Self::UnusedIdent(v)
    }
}

impl ReportSimple for Warnings {
    fn make_simple_report(&self) -> String {
        match self {
            Warnings::UnusedIdent(e) => e.make_simple_report(),
        }
    }
}

impl ReportSimple for LammpsError {
    fn make_simple_report(&self) -> String {
        match self {
            LammpsError::SyntaxError(e) => e.make_simple_report(),
            LammpsError::InvalidStyle(e) => e.make_simple_report(),
            LammpsError::UndefinedIdent(e) => e.make_simple_report(),
            LammpsError::InvalidArguments(e) => e.make_simple_report(),
        }
    }
}

impl From<LammpsError> for lsp_types::Diagnostic {
    fn from(value: LammpsError) -> Self {
        match value {
            LammpsError::SyntaxError(e) => e.into(),
            LammpsError::InvalidStyle(e) => e.into(),
            LammpsError::UndefinedIdent(e) => e.into(),
            LammpsError::InvalidArguments(e) => e.into(),
        }
    }
}

impl From<InvalidStyle> for LammpsError {
    fn from(v: InvalidStyle) -> Self {
        Self::InvalidStyle(v)
    }
}

impl From<UndefinedIdent> for LammpsError {
    fn from(v: UndefinedIdent) -> Self {
        Self::UndefinedIdent(v)
    }
}

impl From<SyntaxError> for LammpsError {
    fn from(error: SyntaxError) -> Self {
        LammpsError::SyntaxError(error)
    }
}

impl From<InvalidArguments> for LammpsError {
    fn from(v: InvalidArguments) -> Self {
        Self::InvalidArguments(v)
    }
}

use crate::diagnostic_report::ReportSimple;
use crate::diagnostics::{Diagnostic, Issue};
use crate::spans::{Point, Span};
use lsp_types::DiagnosticSeverity;
use owo_colors::OwoColorize;
use thiserror::Error;

use crate::fix_styles::FixStyle;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("{}:{}: {err_type} for fix {fix_style}",range.start.row+1,range.start.column+1)]
pub struct InvalidArguments {
    pub(crate) err_type: InvalidArgumentsType,
    pub(crate) range: Span,
    pub(crate) fix_style: FixStyle,
}

impl InvalidArguments {
    pub fn new(err_type: InvalidArgumentsType, range: Span, fix_style: FixStyle) -> Self {
        InvalidArguments {
            err_type,
            range,
            fix_style,
        }
    }

    pub fn range(&self) -> Span {
        self.range
    }
    pub fn start(&self) -> Point {
        self.range.start
    }
    pub fn end(&self) -> Point {
        self.range.end
    }
}

impl Issue for InvalidArguments {
    fn diagnostic(&self) -> crate::diagnostics::Diagnostic {
        Diagnostic {
            name: "Invalid Arguments".to_owned(),
            severity: crate::diagnostics::Severity::Error,
            span: self.range,
            message: format!("for fix `{}`: {}", self.fix_style, self.err_type),
        }
    }
}

impl From<InvalidArguments> for lsp_types::Diagnostic {
    fn from(value: InvalidArguments) -> Self {
        lsp_types::Diagnostic {
            range: value.range.into_lsp_types(),
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

#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum InvalidArgumentsType {
    #[error("Incorrect number of arguments: {provided} provided, expected {expected}")]
    IncorrectNumberArguments { provided: usize, expected: usize },
    #[error(
        "Incorrect keyword arguments for `{kwarg}`. Expected {n_expected} arguments: {expected}. Only {n_provided} provided."
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

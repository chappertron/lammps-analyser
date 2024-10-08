use crate::diagnostics::{Diagnostic, Issue};
use crate::spans::{Point, Span};
use lsp_types::DiagnosticSeverity;
use thiserror::Error;

use crate::styles::FixStyle;

use super::fixes::CommandAndStyle;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
#[error("{err_type} for {style}")]
pub struct InvalidArguments {
    pub(crate) err_type: Box<InvalidArgumentsType>,
    pub(crate) range: Span,
    pub(crate) style: CommandAndStyle,
}

impl InvalidArguments {
    pub fn new(err_type: InvalidArgumentsType, range: Span, style: CommandAndStyle) -> Self {
        InvalidArguments {
            err_type: Box::new(err_type),
            range,
            style,
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
            name: "invalid arguments",
            severity: crate::diagnostics::Severity::Error,
            span: self.range,
            message: format!("for {} : {}", self.style, self.err_type),
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

#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum InvalidArgumentsType {
    #[error("incorrect number of arguments: {provided} provided, expected {expected}")]
    IncorrectNumberArguments { provided: usize, expected: usize },
    #[error(
        "incorrect arguments for keyword `{kwarg}`. expected {n_expected} arguments: {expected}. {n_provided} provided."
    )]
    MissingKwargField {
        kwarg: String,
        expected: String,
        n_expected: usize,
        n_provided: usize,
    },
    #[error("incorrect argument type. expected: {expected}. provided: {provided}.")]
    IncorrectType { expected: String, provided: String },
    #[error("{0}")]
    /// Error for specific fixes
    Custom(String),

    #[error("{} is not a valid keyword for fix {}", kwarg, fix_style)]
    InvalidKeyword { kwarg: String, fix_style: FixStyle },
    #[error("Invalid option `{}` for keyword `{}`. Valid options: [{}].", kwarg, provided,options.join(","))]
    InvalidOption {
        kwarg: &'static str,
        provided: String,
        options: Vec<&'static str>,
    },
    #[error("Invalid style for fix: {0}")]
    InvalidStyle(String),
}

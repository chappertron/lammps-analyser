use crate::diagnostic_report::ReportSimple;
use crate::utils::point_to_position;
use lsp_types::DiagnosticSeverity;
use owo_colors::OwoColorize;
use thiserror::Error;

use tree_sitter::Point;

use crate::fix_styles::FixStyle;

use tree_sitter::Range;

#[derive(Debug, Error, PartialEq, Eq)]
#[error("{}:{}: {err_type} for fix: {fix_style}",range.start_point.row+1,range.start_point.column+1)]
pub struct InvalidArguments {
    pub(crate) err_type: InvalidArgumentsType,
    pub(crate) range: Range,
    pub(crate) fix_style: FixStyle,
}

impl InvalidArguments {
    pub fn new(err_type: InvalidArgumentsType, range: Range, fix_style: FixStyle) -> Self {
        InvalidArguments {
            err_type,
            range,
            fix_style,
        }
    }

    pub fn range(&self) -> Range {
        self.range
    }
    pub fn start(&self) -> Point {
        self.range.start_point
    }
    pub fn end(&self) -> Point {
        self.range.end_point
    }
}

impl From<InvalidArguments> for lsp_types::Diagnostic {
    fn from(value: InvalidArguments) -> Self {
        lsp_types::Diagnostic {
            range: lsp_types::Range {
                start: point_to_position(&(value.start())),
                end: point_to_position(&(value.end())),
            },
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

#[derive(Debug, Error, PartialEq, Eq)]
pub enum InvalidArgumentsType {
    #[error("Incorrect number of arguments: {provided} expected {expected}")]
    IncorrectNumberArguments { provided: usize, expected: usize },
    #[error(
        "Incorrect keyword arguments for `{kwarg}`. 
         Expected {n_expected} arguments: {expected}. 
         Only {n_provided} provided."
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

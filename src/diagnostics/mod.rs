//! Defines the `Issue` trait, the `Diagnostic` struct and the `Severity` Enum
//!
//! These are used for creating a consistent format for found errors to be reported with.

use owo_colors::OwoColorize;
use std::fmt::Display;

use crate::{diagnostic_report::FileNameReport, spans::Span};

/// Indicate that the implement represents an issue that within a LAMMPS script.
///
/// Has a required diagnostic method that gives information about the issue.
pub trait Issue: Sized {
    fn diagnostic(&self) -> Diagnostic;
}

/// A struct to describe a thing as a diagnostic
#[derive(Default, Clone, Eq, PartialEq, Debug)]
pub struct Diagnostic {
    /// Name of the diagnostic.
    pub name: &'static str,
    /// How bad is the diagnostic?
    pub severity: Severity,
    /// The span of the input script responsible for the diagnostic.
    pub span: Span,
    /// Message of the diagnostic.
    pub message: String,
    // Extra information that can be used to understand the issue.
    // TODO: Re add this field once we work out how to add the URI to it..
    // pub information: Vec<Info>,
}

// #[derive(Default, Clone, Eq, PartialEq, Debug)]
// pub struct Info {
//     pub message: String,
//     pub span: Span,
// }

impl FileNameReport for Diagnostic {
    fn make_file_name_report(&self, filename: &str) -> String {
        let start = self.span.start;

        format!(
            "{}: {} in {}:{}:{}",
            self.severity.coloured_display().bold(),
            self.message.bold(),
            filename,
            start.row + 1,
            start.column + 1,
        )
    }
}

impl Display for Diagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let start = self.span.start;
        let end = self.span.end;
        write!(
            f,
            "{}: {}:{}-{}:{},{}",
            self.severity,
            self.message,
            start.row + 1,
            start.column + 1,
            end.row + 1,
            end.column + 1
        )
    }
}

/// The severity of an `Issue` or it's associated `Diagnostic`
///
/// # A Note on `Ord`
/// `Ord` is implemented so Hint < Info < Warning < Error.
/// This means however that the discriminants are different to the `lsp_type::DiagnosticSeverity`
/// values, which go from low to high for most severe to least severe.
#[derive(Default, Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Ord)]
pub enum Severity {
    Hint,
    Info,
    #[default]
    Warning,
    Error,
}

impl From<Severity> for lsp_types::DiagnosticSeverity {
    fn from(value: Severity) -> Self {
        use lsp_types::DiagnosticSeverity as DS;
        match value {
            Severity::Hint => DS::HINT,
            Severity::Info => DS::INFORMATION,
            Severity::Warning => DS::WARNING,
            Severity::Error => DS::ERROR,
        }
    }
}

impl Severity {
    fn coloured_display(&self) -> String {
        let repr = self.to_string().to_lowercase();
        match self {
            Self::Hint => format!("{}", repr.bright_green()),
            Self::Info => format!("{}", repr.bright_blue()),
            Self::Warning => format!("{}", repr.bright_yellow()),
            Self::Error => format!("{}", repr.bright_red()),
        }
    }
}

impl Display for Severity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Hint => write!(f, "HINT"),
            Self::Info => write!(f, "INFO"),
            Self::Warning => write!(f, "WARNING"),
            Self::Error => write!(f, "ERROR"),
        }
    }
}

impl From<Diagnostic> for lsp_types::Diagnostic {
    fn from(value: Diagnostic) -> Self {
        let Diagnostic {
            severity,
            span,
            message,
            ..
        } = value;

        lsp_types::Diagnostic {
            range: span.into_lsp_types(),
            severity: Some(severity.into()),
            source: Some("lammps-analyser".into()),
            message,
            ..Default::default()
        }
    }
}

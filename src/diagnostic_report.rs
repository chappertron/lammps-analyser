use std::fmt::Display;

/// Trait to impl on errors for colourising the output
/// Meant to be for simple output to the commandline
pub trait ReportSimple {
    /// Simply output a string that is printed to the screen.
    fn make_simple_report(&self) -> String;
}

pub trait LspDiagnostic: Display {
    /// Convert the error into an LSP Diagnostic
    /// If the report requires more information than is on the type and the source code,
    /// it may be better to add an inherent method that creates the diagnostic.
    fn make_lsp_report(&self, text: &str) -> lsp_types::Diagnostic;
}

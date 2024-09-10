use std::fmt::Display;

pub trait FileNameReport {
    /// Simply output a string that is printed to the screen.
    fn make_file_name_report(&self, filename: &str) -> String;
}

pub trait LspDiagnostic: Display {
    /// Convert the error into an LSP Diagnostic
    /// If the report requires more information than is on the type and the source code,
    /// it may be better to add an inherent method that creates the diagnostic.
    fn make_lsp_report(&self, text: &str) -> lsp_types::Diagnostic;
}

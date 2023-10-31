use std::fmt::Display;

use ariadne::{Report, Span};
// use std::fmt::Debug;
// use std::hash::Hash;

/// Convert an error into a report that can be outputted into terminal
pub trait ReportDiagnostic<S: Span> {
    fn make_report(&self, source_id: &str) -> Report<'_, S>;
    // where
    //     Id: Debug + Hash + PartialEq + Eq + ToOwned,
    //     <Id as ToOwned>::Owned: From<Id>;
}

// trait LspDiagnostic {
//     fn into_(&self) -> Report;
// }
//
//

/// Trait to impl on errors for colourising the output
pub trait ReportSimple {
    /// Simply output a string that is printed to the screen.
    fn make_simple_report(&self) -> String;
}

pub trait LspDiagnostic: Display {
    /// Convert the error into an LSP Diagnostic
    fn make_lsp_report(&self, text: &str) -> lsp_types::Diagnostic;
}

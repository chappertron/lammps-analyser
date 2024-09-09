use std::{error::Error, fmt::Display};

use crate::{diagnostic_report::ReportSimple, spans::Span};
use thiserror::Error;

/// A simple type that adds an associated `Span` to an error.
/// This information is used to generate better reports, without having to add a span to every
/// possible error type.
///
/// Display is implemented for any E that implements display and transparently passes it through.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Error, Clone, Hash, Copy)]
pub struct SpannedError<E> {
    pub error: E,
    pub span: Span,
}

impl<E> SpannedError<E> {
    pub fn new(error: E, span: impl Into<Span>) -> Self {
        SpannedError {
            error,
            span: span.into(),
        }
    }
}

impl<V> SpannedError<V> {
    /// Convert the inner error, using the appropriate `From` impl.
    pub fn convert<U>(value: SpannedError<U>) -> Self
    where
        V: From<U>,
    {
        Self {
            error: value.error.into(),
            span: value.span,
        }
    }
}

/// Add a span to the error type.
pub(crate) trait WithSpan {
    type Output;
    type Error;
    fn with_span(self, span: Span) -> Result<Self::Output, SpannedError<Self::Error>>;
}

impl<T, E: Error> WithSpan for Result<T, E> {
    type Output = T;
    type Error = E;
    fn with_span(self, span: Span) -> Result<Self::Output, SpannedError<Self::Error>> {
        match self {
            Ok(x) => Ok(x),
            Err(error) => Err(SpannedError { error, span }),
        }
    }
}

impl<E> Display for SpannedError<E>
where
    E: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl<E> ReportSimple for SpannedError<E>
where
    E: ReportSimple,
{
    fn make_simple_report(&self) -> String {
        let start = self.span.start;
        format!(
            "{}:{}: {}",
            start.row,
            start.column,
            self.error.make_simple_report()
        )
    }
}

// TODO: Implement the report trait??

use crate::spans::Span;
use thiserror::Error;

/// A simple type that adds an associated `Span` to an error.
/// This information is used to generate better reports, without having to add a span to every
/// possible error type.
#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord, Error)]
pub struct SpannedError<E> {
    pub error: E,
    pub span: Span,
}

// TODO: Implement the report trait??

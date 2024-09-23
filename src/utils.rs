use crate::spans::Point;

use crate::identifinder::{IdentiFinder, NameAndType};

pub(crate) mod tree_sitter_helpers;

/// Find the symbol under the given point in the file
pub(crate) fn get_symbol_at_point<'a>(
    identifinder: &'a IdentiFinder,
    point: &Point,
) -> Option<&'a NameAndType> {
    let symbols = identifinder.symbols();

    // Might've created something a little slow here...
    for (k, v) in symbols {
        for r in v.refs().iter() {
            if *point >= r.start() && *point <= r.end() {
                return Some(k);
            }
        }
    }

    None
}

pub(crate) mod into_error {

    /// A trait designed for allowing easier conversions from `Option<T>` to `Result<T, E>`.
    ///
    /// This needed for foreign types T to circumvent the orphan rules, that prevent using the `From` trait.
    /// This can only be defined Once per type, but then `From<Error>` can be implemented for other
    /// error types with the question mark operator.
    pub(crate) trait IntoError<T> {
        type Error;
        // TODO: allow for custom messages
        fn into_err(self) -> Result<T, Self::Error>;
    }
}

pub(crate) mod parsing {
    use tree_sitter::Parser;

    /// Create an instance of the tree-sitter-lammps parser
    pub fn setup_parser() -> Parser {
        let mut parser = Parser::new();
        parser.set_language(tree_sitter_lammps::language()).unwrap();
        parser
    }
}

#[cfg(test)]
pub(crate) mod testing {
    use tree_sitter::Parser;

    pub fn setup_parser() -> Parser {
        super::parsing::setup_parser()
    }
    pub(crate) fn parse(source_bytes: impl AsRef<[u8]>) -> tree_sitter::Tree {
        let source_bytes = source_bytes.as_ref();
        let mut parser = setup_parser();
        parser.parse(source_bytes, None).unwrap()
    }
}

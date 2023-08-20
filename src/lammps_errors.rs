use crate::error_finder::SyntaxError;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum LammpsError {
    #[error("{0}")]
    SyntaxError(SyntaxError),
}

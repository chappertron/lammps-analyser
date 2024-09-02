//! Checks the validity of fix arguments

use thiserror::Error;

use crate::{
    ast::{Ast, CommandType, NamedCommand},
    commands::CommandName,
    diagnostics::{self, Diagnostic, Issue},
    spans::Span,
};

use self::{computes::check_compute, fixes::check_fix, invalid_arguments::InvalidArguments};

pub mod invalid_arguments;

pub mod computes;
pub mod fixes;

pub mod utils;

#[derive(Debug, Clone, Eq, PartialEq, Error)]
pub enum InvalidCommand {
    #[error("Unknown command: {0}")]
    UnknownCommand(String, Span),
    #[error("{0}")]
    InvalidArguments(InvalidArguments),
}

impl From<InvalidArguments> for InvalidCommand {
    fn from(v: InvalidArguments) -> Self {
        Self::InvalidArguments(v)
    }
}

impl Issue for InvalidCommand {
    fn diagnostic(&self) -> diagnostics::Diagnostic {
        let span = match self {
            Self::UnknownCommand(_, span) => *span,
            Self::InvalidArguments(invalid_args) => invalid_args.range,
        };
        Diagnostic {
            name: "Invalid Command".to_owned(),
            message: self.to_string(),
            span,
            severity: diagnostics::Severity::Error,
        }
    }
}

impl Ast {
    pub fn check_commands(&self) -> impl Iterator<Item = InvalidCommand> + '_ {
        self.commands
            .iter()
            .filter_map(|command| match &command.command_type {
                CommandType::NamedCommand(NamedCommand::Fix(fix_def)) => {
                    // Check if the fixes are ok.
                    check_fix(fix_def).err().map(InvalidCommand::from)
                }
                CommandType::GenericCommand(command) => {
                    let command_name = CommandName::from(command.name.as_str());

                    if let CommandName::InvalidCommand(name) = command_name {
                        Some(InvalidCommand::UnknownCommand(
                            name.clone(),
                            Span {
                                start: command.start,
                                end: command.end,
                            },
                        ))
                    } else {
                        None
                    }
                }
                // Not sure about other named commands yet...
                CommandType::NamedCommand(NamedCommand::Compute(compute)) => {
                    check_compute(compute).err().map(InvalidCommand::from)
                }
                _ => todo!(),
            })
    }
}

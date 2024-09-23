use itertools::Itertools;

use crate::spans::{Point, Span};

use crate::ast;

impl ast::Ast {
    /// Find the command corresponding to a given point
    pub fn find_point(&self, point: &Point) -> Option<&ast::Command> {
        self.commands
            .iter()
            .find(|cmd| cmd.span().start <= *point && cmd.span().end >= *point)
    }

    /// An iterator over references to run commands
    pub fn find_run_commands(&self) -> impl Iterator<Item = &ast::Command> {
        // TODO: also include `minimize` and `rerun` commands
        self.commands.iter().filter(|cmd| {
            if let ast::Command::GenericCommand(c) = &cmd {
                c.name.contents == "run"
            } else {
                false
            }
        })
    }

    // TODO: Is this lifetime bound correct?
    /// Find the spans of commands before a each run command.
    ///
    /// These spans start at the end of the previous run command and end at the start of the
    /// current run command
    ///
    /// The first span starts at line zero column zero
    pub fn find_run_blocks(&self) -> impl Iterator<Item = Span> + '_ {
        let run_definitions = self.find_run_commands();
        // spans of the run command themselves
        let run_spans = run_definitions.map(|cmd| cmd.span());

        // The spans of the blocks of commands defined before run commands.
        std::iter::once(Span {
            start: Default::default(),
            end: Default::default(),
        })
        .chain(run_spans)
        .tuple_windows()
        .map(|(cmd1, cmd2)| Span {
            start: cmd1.end,
            end: cmd2.start,
        })
    }
}

#[cfg(test)]
mod test {
    use crate::input_script::InputScript;
    use pretty_assertions::assert_eq;

    use super::*;
    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_find_run_commands() {
        let file = include_str!("../../example_input_scripts/in.nemd");

        let input_script = InputScript::new(file).unwrap();

        let run_commands: Vec<_> = input_script.ast.find_run_commands().cloned().collect();

        assert_eq!(run_commands.len(), 2);
        assert_eq!(run_commands[0].span(), Span::from((131, 0)..(131, 17)));
        assert_eq!(run_commands[1].span(), Span::from((145, 0)..(145, 6)));
    }
}

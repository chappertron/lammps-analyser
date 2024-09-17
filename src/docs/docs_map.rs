use std::collections::HashMap;

use once_cell::sync::Lazy;

use crate::{commands::CommandName, compute_styles::ComputeStyle, fix_styles::FixStyle};
/// Read the `index_map' file generated by the docs cleaning up code.
/// This maps between commands and the documentation files

/// Static map between the commands and their files.
/// TODO: Make pub(crate) again
pub static DOCS_MAP: Lazy<DocsMap> =
    Lazy::new(|| DocsMap::new_from_str(include_str!["../../docs_extract/index_map.txt"]));

/// A mapping between a command name and a doc path.
pub struct DocsMap {
    fixes: HashMap<FixStyle, String>,
    computes: HashMap<ComputeStyle, String>,
    commands: HashMap<CommandName, String>,
}

impl DocsMap {
    /// Read an 'index file' that maps commands to their
    ///
    /// # Panics
    ///
    /// Currently panics if the 'map' file does not have comma separated pairs of
    /// index names and file prefixes.
    pub fn new_from_str(map_contents: &str) -> Self {
        let mut fixes = HashMap::new();
        let mut computes = HashMap::new();
        let mut others = HashMap::new();

        map_contents.lines().for_each(|line| {
            // TODO: Remove this panic.
            let (name, doc_file) = line.split_once(',').expect("Expected key value pair");

            if let Some((command_type, style)) = name.split_once(' ') {
                match command_type {
                    "fix" => {
                        fixes.insert(style.into(), doc_file.to_owned());
                    }
                    "compute" => {
                        computes.insert(style.into(), doc_file.to_owned());
                    }
                    _ => {
                        others.insert(style.into(), doc_file.to_owned());
                    }
                }
            } else {
                // Name can't be split, so obviously has no style
                others.insert(name.into(), doc_file.to_owned());
            }
        });

        DocsMap {
            fixes,
            computes,
            commands: others,
        }
    }

    /// Map between fix styles and their appropriate file paths
    pub fn fixes(&self) -> &HashMap<FixStyle, String> {
        &self.fixes
    }

    /// Map between compute styles and their appropriate file paths
    pub fn computes(&self) -> &HashMap<ComputeStyle, String> {
        &self.computes
    }

    /// Map between commands and their appropriate file paths
    pub fn commands(&self) -> &HashMap<CommandName, String> {
        &self.commands
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        commands::CommandName, compute_styles::ComputeStyle, docs::docs_map::DOCS_MAP,
        fix_styles::FixStyle,
    };

    use super::DocsMap;

    #[test]
    fn index_map() {
        let file = include_str!("../../docs_extract/index_map.txt");

        let map = DocsMap::new_from_str(file);

        // Obvious example
        assert_eq!(map.fixes()[&FixStyle::Nve], "fix_nve".to_owned());

        // Non-obvious example
        assert_eq!(map.fixes()[&FixStyle::Nvt], "fix_nh".to_owned());

        // Obvious example
        assert_eq!(
            map.computes()[&ComputeStyle::StressAtom],
            "compute_stress_atom".to_owned()
        );

        // Obvious example
        assert_eq!(map.commands()[&CommandName::Newton], "newton".to_owned());
    }

    #[test]
    fn index_map_static() {
        let map = &DOCS_MAP;

        // Obvious example
        assert_eq!(map.fixes()[&FixStyle::Nve], "fix_nve".to_owned());

        // Non-obvious example
        assert_eq!(map.fixes()[&FixStyle::Nvt], "fix_nh".to_owned());

        // Obvious example
        assert_eq!(
            map.computes()[&ComputeStyle::StressAtom],
            "compute_stress_atom".to_owned()
        );

        // Obvious example
        assert_eq!(map.commands()[&CommandName::Newton], "newton".to_owned());
    }
}

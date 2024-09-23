# LAMMPS-Analyser Changelog

## Initial Release 0.1.0

### Features

- A cli `lammps-analyser` that checks input scripts
- A language server `lmp-lsp` that performs these checks within the text editor.
- Supported checks:
  - Warnings about unused variables
  - Errors if variables/fixes/computes are used but not defined.
  - Commands have valid names
  - That fix/compute/pair styles are valid.
  - The number of arguments for fix and compute commands.
  - Fully checks arguments for the npt/nph/npt fixes.
  - Warns if computes fixes or variables have been defined multiple times before a `run` command. 
- The LSP provides:
  - Go-to definition and references for variables, computes and fixes.
  - Hover documentation for commands, fixes and computes.
  - Checks your script as you type and reports errors and warnings 
- Based of LAMMPS 27June2024 Feature Release

### Internal Features

- Scripts to convert LAMMPS documentation files to Markdown and extract `fix`,`pair` and `compute` styles
- A parser that converts a `tree-sitter` tree of the script into an abstract syntax tree.


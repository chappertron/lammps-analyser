# Release Checklist

## lammps-analyser

- [x] Bump version (to 0.1.0)
- [x] Change dependency on `tree-sitter-lammps` to be via git.
- [x] Change dependency on `tree-sitter-lammps` to be via crates.
- [x] Add `README.md`.
- [ ] Add pictures to finish the read-me
- [x] Add licence
- [ ] Publish to crates.io
- [x] Add CI for testing
- [x] Add CI for release binaries (`cargo-dist`)
- [ ] Create a release tag `gh release`
- [ ] Test release CI with a pre-release
- [x] Fix warnings in compiling
- [x] Run Clippy
- [x] Run typos
- [x] Check that all docs are committed and builds in CI.
- [x] Rename `lsp` to `lmp-lsp`
- [x] Create a CHANGELOG.md
- [x] Allow for quoted expressions in variable commands

## tree-sitter lammps

- [x] Publish to crates.io
- [x] Fix warnings in compiling the grammar in RUST.

# Release Checklist

## lammps-analyser

- [ ] Bump version (to 0.1.3-beta if possible)
- [x] Change dependency on treesitter-lammps to be via git.
- [ ] Change dependency on treesitter-lammps to be via crates.
- [ ] Add README.md
- [ ] Add licence
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
- [ ] Create a CHANGELOG.md
- [x] Allow for quoted expressions in variable commands

## tree-sitter lammps

- [ ] Publish to crates.io
- [ ] Investigate what should be committed, especially re. WASM.
- [x] Fix warnings in compiling the grammar in RUST.

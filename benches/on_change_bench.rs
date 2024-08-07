use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lammps_analyser::lsp;
use tower_lsp::Client;

pub fn lsp_on_change(c: &mut Criterion) {
    todo!()
}

criterion_group!(benches, lsp_on_change);
criterion_main!(benches);

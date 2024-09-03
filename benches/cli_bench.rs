use criterion::{black_box, criterion_group, criterion_main, Criterion};
use lammps_analyser::input_script::InputScript;

pub fn read_script_w_errors(c: &mut Criterion) {
    // setup
    let mut parser = tree_sitter::Parser::new();
    parser
        .set_language(tree_sitter_lammps::language())
        .expect("Failed to setup tree-sitter-lammps");

    let source = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/example_input_scripts/in.nemd"
    ));

    c.bench_function("read in.nemd", |b| {
        b.iter(|| {
            black_box(InputScript::new(source, &mut parser).expect("failed to catch errors."));
        })
    });
}

criterion_group!(benches, read_script_w_errors);
criterion_main!(benches);

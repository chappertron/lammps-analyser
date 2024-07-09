use std::{
    env,
    error::Error,
    fs::{self, File},
    io::Write,
    path::Path,
};

// Relative to crate root
const SOURCE_DIR: &str = "lammps_docs_md/";

// TODO: Add rerun checking logic
fn main() -> Result<(), Box<dyn Error>> {
    let out_dir = env::var("OUT_DIR")?;
    let dest_path = Path::new(&out_dir).join("all_the_docs.rs");
    let crate_dir = env::var("CARGO_MANIFEST_DIR")?;
    let crate_root = Path::new(&crate_dir);
    let mut all_the_docs = File::create(&dest_path)?;

    dbg!(out_dir);
    writeln!(&mut all_the_docs, r##"["##,)?;

    for f in fs::read_dir(crate_root.join(SOURCE_DIR))? {
        let f = f?;

        if !f.file_type()?.is_file() {
            continue;
        }
        let name = f
            .file_name()
            .into_string()
            // TODO: Proper error handling
            .expect("Filename is invalid utf-8");

        writeln!(
            &mut all_the_docs,
            r##"("{name}", include_str!(r#"{path}"#)),"##,
            path = f.path().display(),
        )?;
    }

    writeln!(&mut all_the_docs, r##"]"##,)?;

    Ok(())
}

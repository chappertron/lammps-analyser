use std::{
    env,
    error::Error,
    fmt::Display,
    fs::{self, File},
    io::Write,
    path::Path,
};

// Relative to crate root
const SOURCE_DIR: &str = "lammps_docs_md/";

fn main() -> Result<(), Box<dyn Error>> {
    let out_dir = env::var("OUT_DIR")?;
    let dest_path = Path::new(&out_dir).join("all_the_docs.rs");
    let crate_dir = env::var("CARGO_MANIFEST_DIR")?;
    let crate_root = Path::new(&crate_dir);
    let mut all_the_docs = File::create(dest_path)?;

    dbg!(out_dir);
    writeln!(&mut all_the_docs, r##"["##,)?;

    for f in fs::read_dir(crate_root.join(SOURCE_DIR))? {
        let f = f?;

        if !f.file_type()?.is_file() {
            continue;
        }
        let name = f.file_name().into_string().map_err(|_| Utf8Error)?;

        writeln!(
            &mut all_the_docs,
            r##"("{name}", include_str!(r#"{path}"#)),"##,
            path = f.path().display(),
        )?;
    }

    writeln!(&mut all_the_docs, r##"]"##,)?;

    // Tell cargo to only rebuild if the docs files have changed.
    println!("cargo::rerun-if-changed={SOURCE_DIR}");

    Ok(())
}

#[derive(Debug)]
struct Utf8Error;
impl Display for Utf8Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "OS String had invalid UTF-8.")
    }
}
impl Error for Utf8Error {}

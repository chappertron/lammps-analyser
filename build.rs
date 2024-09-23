use std::{
    env,
    error::Error,
    fmt::Display,
    fs::{self, create_dir_all, read_to_string, File},
    io::Write,
    path::{Path, PathBuf},
    str::FromStr,
};

// Relative to crate root
const SOURCE_DIR: &str = "lammps_docs_md/";

fn main() -> BResult<()> {
    let out_dir = env::var("OUT_DIR")?;
    let dest_path = Path::new(&out_dir).join("all_the_docs.rs");
    let crate_dir = env::var("CARGO_MANIFEST_DIR")?;
    let crate_root = Path::new(&crate_dir);
    let all_the_docs = File::create(dest_path)?;

    generate_docs_map(all_the_docs, crate_root)?;

    generate_styles(StyleKind::PairStyle)?;

    // Tell cargo to only rebuild if the docs files have changed.
    println!("cargo::rerun-if-changed={SOURCE_DIR}");

    Ok(())
}

type BResult<O> = Result<O, Box<dyn Error>>;

/// Generate the mapping between the index key and the documentation files.
/// These are included in a 'map'
fn generate_docs_map(mut out_file: File, crate_root: &Path) -> BResult<()> {
    writeln!(&mut out_file, r##"["##,)?;

    for f in fs::read_dir(crate_root.join(SOURCE_DIR))? {
        let f = f?;

        if !f.file_type()?.is_file() {
            continue;
        }
        let name = f.file_name().into_string().map_err(|_| Utf8Error)?;

        writeln!(
            &mut out_file,
            r##"("{name}", include_str!(r#"{path}"#)),"##,
            path = f.path().display(),
        )?;
    }

    writeln!(&mut out_file, r##"]"##,)?;
    Ok(())
}

enum StyleKind {
    PairStyle,
}

struct StyleInfo {
    /// The path of the source file.
    file_name: &'static str,
    /// The name of the generated enum.
    enum_name: &'static str,
    /// Where the generated enum will be stored.
    out_path: &'static str,
}

impl StyleKind {
    fn file_name(&self) -> &str {
        self.info().file_name
    }

    fn enum_name(&self) -> &str {
        self.info().enum_name
    }

    fn out_path(&self) -> &str {
        self.info().out_path
    }

    fn info(&self) -> StyleInfo {
        match self {
            Self::PairStyle => StyleInfo {
                file_name: "pair_styles.txt",
                enum_name: "PairStyle",
                out_path: "pair_styles.rs",
            },
        }
    }
}

fn generate_styles(style_kind: StyleKind) -> BResult<()> {
    let fname = style_kind.file_name();
    let mut path = PathBuf::from_str(&std::env::var("CARGO_MANIFEST_DIR")?)?;
    path.push("docs_extract");
    path.push(fname);

    let styles = read_to_string(path)?;

    let output = Path::new(&env::var("OUT_DIR")?)
        .join("styles")
        .join(style_kind.out_path());

    create_dir_all(output.parent().expect("Should have parent dir"))?;

    let mut output = File::create(output)?;

    // Open the macro
    writeln!(output, "derive_styles! [",)?;
    writeln!(output, "{} =>", style_kind.enum_name())?;

    // Write the enum variants.
    for style in styles.lines() {
        let variant = camelify_lammps(style);
        writeln!(output, r##"({variant},"{style}"),"##)?;
    }

    // Close the macro
    writeln!(output, "];",)?;

    Ok(())
}

/// Convert from 'lammps/case' to CamelCase
///
/// Words that are snake case also get further split
fn camelify_lammps(src: &str) -> String {
    src.split('/').map(camelify_snake).collect()
}

fn camelify_snake(src: &str) -> String {
    src.split('_').map(uppercase_first_letter).collect()
}

fn uppercase_first_letter(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}

#[derive(Debug)]
struct Utf8Error;
impl Display for Utf8Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "OS String had invalid UTF-8.")
    }
}
impl Error for Utf8Error {}

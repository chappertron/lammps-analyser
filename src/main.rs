/// TODO:
/// - [x] Scan all identifiers in one pass?
/// - [ ] Add incrental parsing
/// - [x] Make the output look a bit prettier
/// - [ ] Add file snippets to show where errors are
/// - [ ] Add a test suite
/// - [ ] Create an issue abstraction to handle all errors and warnings
use anyhow::Result;
use clap::Parser as ClapParser;
use lammps_analyser::{
    compute_styles::ComputeStyle, error_finder::ErrorFinder, fix_styles::FixStyle,
    identifinder::IdentiFinder,
};
use owo_colors::OwoColorize;
use std::{
    fmt::Display,
    fs::File,
    io::{BufReader, Read},
};
use thiserror::Error;
use tree_sitter::{Parser, Point, Query, QueryCursor, Tree};

#[derive(Debug, clap::Parser)]

struct Cli {
    source: String,
    /// Output the parsed tree to a dot file
    #[clap(long)]
    output_tree: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let file = File::open(&cli.source)?;

    let mut file = BufReader::new(file);

    let mut source_code = String::new();

    file.read_to_string(&mut source_code)?;

    let source_bytes = source_code.as_bytes();

    let mut parser = Parser::new();

    parser
        .set_language(tree_sitter_lammps::language())
        .expect("Could not load language");

    let tree = parser.parse(source_bytes, None).unwrap();

    // dbg!(&tree.root_node().to_sexp());

    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        tree.print_dot_graph(&dot_file);
    }

    let mut identifinder = IdentiFinder::new(&tree, source_bytes)?;

    let undefined_fixes = match identifinder.check_symbols() {
        Ok(()) => vec![],
        Err(v) => v,
    };

    let mut error_finder = ErrorFinder::new()?;
    _ = error_finder.find_syntax_errors(&tree, source_bytes)?;
    error_finder.find_missing_nodes(&tree)?;
    let syntax_errors = error_finder.syntax_errors();

    for err in syntax_errors {
        // ruff_test.py:3:1: F821 Undefined name `hello`
        println!("{}:{}", cli.source.bold(), err);
    }

    // this doesn't work. Need to compare the names!

    for err in &undefined_fixes {
        // ruff_test.py:3:1: F821 Undefined name `hello`
        // println!("{}",std::str::from_utf8(source_code[ident.start_byte..ident.end_byte])?.underline());
        println!("{}:{}", cli.source.bold(), err);
    }

    let invalid_styles = check_styles(&tree, source_bytes)?;
    for err in &invalid_styles {
        println!("{}:{}", cli.source.bold(), err);
    }
    // TODO Check if any warnings or errors are found!!!

    if !syntax_errors.is_empty() || !undefined_fixes.is_empty() || !invalid_styles.is_empty() {
        let n_errors = syntax_errors.len() + undefined_fixes.len() + invalid_styles.len();
        println!(
            "{}: {} error{} found 😞",
            cli.source.bold(),
            n_errors.bright_red(),
            if n_errors != 1 { "s" } else { "" },
        );
        Ok(())
    } else {
        println!("All Good 😊");
        Ok(())
    }
}

#[derive(Debug, Error)]
#[error("{}:{}: {} {} `{}`",start.row+1,start.column+1,"Invalid".bright_red(),style_type.bright_red(),name)]
pub struct InvalidStyle {
    pub start: Point,
    pub end: Point,
    pub name: String,
    pub style_type: StyleType,
}

#[derive(Debug, PartialEq, Eq)]
pub enum StyleType {
    Fix,
    Compute,
    Pair,
    Kspace,
    Minimize,
}

impl Display for StyleType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fix => write!(f, "fix style"),
            Self::Compute => write!(f, "compute style"),
            Self::Pair => write!(f, "pair style"),
            Self::Kspace => write!(f, "kspace style"),
            Self::Minimize => write!(f, "minimize style"),
        }
    }
}

fn check_styles(tree: &Tree, text: &[u8]) -> Result<Vec<InvalidStyle>> {
    let query = Query::new(
        tree.language(),
        "(fix (fix_style) @definition.fix) (compute (compute_style) @definition.compute) ",
    )?;

    let mut query_cursor = QueryCursor::new();

    let matches = query_cursor.matches(&query, tree.root_node(), text);

    Ok(matches
        .into_iter()
        .filter_map(|mat| {
            let style = mat.captures[0].node.utf8_text(text).ok()?;
            let style_type = match mat.captures[0].node.kind() {
                "fix_style" => StyleType::Fix,
                "compute_style" => StyleType::Compute,
                _ => unreachable!(),
            };

            if let (StyleType::Fix, FixStyle::InvalidFixStyle) = (&style_type, style.into()) {
                Some(InvalidStyle {
                    start: mat.captures[0].node.start_position(),
                    end: mat.captures[0].node.end_position(),
                    name: style.to_string(),
                    style_type: StyleType::Fix,
                })
            } else if let (StyleType::Compute, ComputeStyle::InvalidComputeStyle) =
                (&style_type, style.into())
            {
                Some(InvalidStyle {
                    start: mat.captures[0].node.start_position(),
                    end: mat.captures[0].node.end_position(),
                    name: style.to_string(),
                    style_type: StyleType::Compute,
                })
            } else {
                None
            }
        })
        .collect())
}
// fn find_undefined_fix_names() {}

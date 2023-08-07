/// TODO:
/// - [ ] Scan all identifiers in one pass?
/// - [ ] Add incrental parsing
/// - [ ] Make the output look a bit prettier
/// - [ ] Add file snippets to show where errors are
/// - [ ] Add a test suite
use anyhow::{Ok, Result};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, Read},
};
use tree_sitter::{Parser, Query, QueryCursor, Tree};
use lammps_analyser::identifinder::{Ident};
use clap::Parser as ClapParser;

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

    let source_code = source_code.as_bytes();

    let mut parser = Parser::new();

    parser
        .set_language(tree_sitter_lammps::language())
        .expect("Could not load language");

    let tree = parser.parse(source_code, None).unwrap();

    // dbg!(&tree.root_node().to_sexp());

    if cli.output_tree {
        let dot_file = File::create("tree.dot")?;
        tree.print_dot_graph(&dot_file);
    }

    let used_fixes = get_fix_names(&tree, source_code)?;
    let def_fixes = get_fix_defs(&tree, source_code)?;

    // dbg!(used_fixes
    //     .iter()
    //     .map(|node| (node.utf8_text(source_code).unwrap(), node.id()))
    //     .collect::<Vec<_>>());

    // this doesn't work. Need to compare the names!
    for ident in used_fixes.difference(&def_fixes) {
        // ruff_test.py:3:1: F821 Undefined name `hello`
        println!(
            "{}:{}:{}: Undefined fix `{}`",
            cli.source,
            ident.start.row + 1,
            ident.start.column + 1,
            ident.name,
        );
    }

    Ok(())
}

fn get_fix_names(tree: &Tree, text: &[u8]) -> Result<HashSet<Ident>> {
    let query = Query::new(tree.language(), "(fix_id) @reference.fix")?;

    let mut query_cursor = QueryCursor::new();

    let matches = query_cursor.matches(&query, tree.root_node(), text);

    // for capture in captures {
    // dbg!(capture.0.captures[0].node.utf8_text(text)?);
    // }

    let fix_ids = matches
        // .map(|mat| mat.captures[0].node.utf8_text(text))
        .map(|mat| Ident::new(&mat.captures[0].node, text))
        // .collect::<HashSet<_>>();
        .collect::<Result<HashSet<_>>>()?;

    Ok(fix_ids)
}
/// List of fixs names that have beend defined in the current scope
fn get_fix_defs(tree: &Tree, text: &[u8]) -> Result<HashSet<Ident>> {
    let query = Query::new(tree.language(), "(fix (fix_id ) @definition.fix) ")?;
    

    let mut query_cursor = QueryCursor::new();

    let matches = query_cursor.matches(&query, tree.root_node(), text);
    



    let fix_ids = matches
        // .map(|mat| mat.captures[0].node.utf8_text(text))
        .map(|mat| Ident::new(&dbg![mat.captures[0]].node, text))
        // .collect::<HashSet<_>>();
        .collect::<Result<HashSet<_>>>()?;

    Ok(fix_ids)
}

// fn find_undefined_fix_names() {}


//! Convert the treesiterr trees into an AST.
//! TODO Work out why the generic command nodes are causing errors:

use tree_sitter::{Node, Tree};

use crate::identifinder::Ident;

pub fn ts_to_ast(tree: &Tree, text: &[u8]) -> Ast {
    let mut cursor = tree.walk();

    // println!("{}", tree.root_node().to_sexp());
    //while cursor.goto_first_child() {
    cursor.goto_first_child();
    println!("{}", cursor.node().to_sexp());
    loop {
        // println!("Current position {:?}", cursor.node().start_position());
        // dbg!(cursor.node().children(&mut cursor).collect::<Vec<_>>());
        // println!("{}", cursor.node().to_sexp());
        if cursor.goto_first_child() && cursor.node().kind() != "comment" {
            println!("{}", cursor.node().to_sexp());
            match NamedCommand::try_from(cursor.node().kind()) {
                Ok(cmd) => {
                    println!("Found named command: {:?}", cmd);
                }
                Err(_) => {
                    println!(
                        "Found generic command: {:?}",
                        GenericCommand::from_node(&cursor.node(), &mut cursor, text)
                    );
                }
            }
            cursor.goto_parent();
        }

        if !cursor.goto_next_sibling() {
            break;
        }
    }

    //}

    todo!()
}

pub struct Ast {
    pub commands: Vec<Command>,
}

pub enum Command {
    GenericCommand(GenericCommand),
    NamedCommand(NamedCommand),
}

#[derive(Debug)]
pub struct GenericCommand {
    pub name: String,
    pub args: Vec<Argument>,
}

impl GenericCommand {
    fn from_node(
        node: &Node,
        cursor: &mut tree_sitter::TreeCursor,
        text: &[u8],
    ) -> Result<Self, String> {
        // let kind = node.kind().to_string();
        let mut args = vec![];
        println!("Debug GenericCommand");

        dbg!(node.to_sexp());
        dbg!(cursor.goto_first_child());
        dbg!(cursor.goto_first_child());

        // TODO use a field in the TS grammar
        let name = cursor.node().utf8_text(text).unwrap().to_string();
        dbg!(&name);
        while dbg!(cursor.goto_next_sibling()) {
            dbg!(node.to_sexp());
            args.push(Argument::from_node(node, cursor, text)?);
        }

        cursor.goto_parent();
        Ok(GenericCommand { name, args })
    }
}

/// Acceptable argument types for LAMMPS commands
#[derive(Debug)]
pub enum Argument {
    /// Expression from LAMMPS
    Expression,
    String,
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    /// TODO make this hold and `Ident` struct, from `crate::identifinder` module
    UnderscoreIdent(Ident),
    /// Variables within curly braces
    /// TODO make this hold and `Ident` struct, from `crate::identifinder` module
    CurlyVar(Ident),
}

impl Argument {
    fn from_node(
        node: &Node,
        _cursor: &mut tree_sitter::TreeCursor,
        text: &[u8],
    ) -> Result<Self, String> {
        // TODO make these variants more complete.
        match node.kind() {
            "expression" => Ok(Self::Expression),
            "string" => Ok(Self::String),
            "group" => Ok(Self::Group),
            "underscore_ident" => Ok(Self::UnderscoreIdent(
                Ident::new(node, text).map_err(|x| format!("{}", x))?,
            )),
            "curly_var" => Ok(Self::CurlyVar(
                Ident::new(node, text).map_err(|x| format!("{}", x))?,
            )),

            x => Err(format!("Unknown argument type: {}", x)),
        }
    }
}

/// Commands that have a special form in the tree sitter grammar
#[derive(Debug)]
pub enum NamedCommand {
    Fix,
    Compute,
    Style,
    Modify,
    AtomStyle,
    Boundary,
    VariableDef,
    ThermoStyle,
    Thermo,
    Units,
    Run,
}

impl TryFrom<&str> for NamedCommand {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "fix" => Ok(Self::Fix),
            "compute" => Ok(Self::Compute),
            "style" => Ok(Self::Style),
            "modify" => Ok(Self::Modify),
            "atom_style" => Ok(Self::AtomStyle),
            "boundary" => Ok(Self::Boundary),
            "variable" => Ok(Self::VariableDef),
            "thermo_style" => Ok(Self::ThermoStyle),
            "thermo" => Ok(Self::Thermo),
            "units" => Ok(Self::Units),
            "run" => Ok(Self::Run),
            _ => Err("Not a named command".into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use tree_sitter::Parser;

    use super::ts_to_ast;

    #[test]
    fn test_ast() {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        let source_bytes = include_bytes!("../../fix.lmp");
        let tree = parser.parse(source_bytes, None).unwrap();

        let ast = ts_to_ast(&tree, source_bytes);
        unimplemented!()
    }
}

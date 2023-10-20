//! Convert the treesiterr trees into an AST.

pub mod expressions;
use tree_sitter::{Node, Point, Tree};

use crate::identifinder::Ident;

pub fn ts_to_ast(tree: &Tree, text: &[u8]) -> Ast {
    let mut cursor = tree.walk();

    let mut commands = vec![];
    // println!("{}", tree.root_node().to_sexp());
    //while cursor.goto_first_child() {
    cursor.goto_first_child();
    println!("{}", cursor.node().to_sexp());
    loop {
        // println!("Current position {:?}", cursor.node().start_position());
        // dbg!(cursor.node().children(&mut cursor).collect::<Vec<_>>());
        // println!("{}", cursor.node().to_sexp());

        // Advance cursor and skip if a comment
        if cursor.goto_first_child() && cursor.node().kind() != "comment" {
            println!("{}", cursor.node().to_sexp());
            // match NamedCommand::try_from(cursor.node().kind()) {
            //     Ok(cmd) => {
            //         println!("Found named command: {:?}", cmd);
            //     }
            //     Err(_) => {
            //         println!(
            //             "Found generic command: {:?}",
            //             GenericCommand::from_node(&cursor.node(), &mut cursor, text)
            //         );
            //     }
            // }
            //

            let cmd = if let Ok(cmd) = NamedCommand::try_from(cursor.node().kind()) {
                // TODO add arguments
                Command::NamedCommand(cmd)
            } else {
                Command::GenericCommand(
                    GenericCommand::from_node(&cursor.node(), &mut cursor, text).unwrap(),
                )
            };

            commands.push(cmd);

            cursor.goto_parent();
        }

        // If no more commands, break!
        if !cursor.goto_next_sibling() {
            break;
        }
    }

    Ast { commands }
}

#[derive(Debug)]
/// A list of commands
/// Perhaps not truly an AST
pub struct Ast {
    pub commands: Vec<Command>,
}

#[derive(Debug)]
pub enum Command {
    GenericCommand(GenericCommand),
    NamedCommand(NamedCommand),
}

#[derive(Debug)]
pub struct GenericCommand {
    pub name: String,
    pub args: Vec<Argument>,
    pub start: Point,
    pub end: Point,
    pub start_byte: usize,
    pub end_byte: usize,
}

impl GenericCommand {
    fn from_node(
        node: &Node,
        cursor: &mut tree_sitter::TreeCursor,
        text: &[u8],
    ) -> Result<Self, String> {
        // let kind = node.kind().to_string();
        let start = node.start_position();
        let end = node.end_position();
        let start_byte = node.start_byte();
        let end_byte = node.end_byte();

        let mut args = vec![];
        println!("Debug GenericCommand");

        assert!(cursor.node() == *node);

        cursor.goto_first_child();

        // TODO use a field in the TS grammar
        let name = cursor.node().utf8_text(text).unwrap().to_string();
        // dbg!(&name);
        while cursor.goto_next_sibling() {
            // dbg!(cursor.node().to_sexp());
            for node in cursor.node().children(cursor) {
                args.push(Argument::from_node(&node, text)?);
            }
        }

        cursor.goto_parent();
        Ok(GenericCommand {
            name,
            args,
            start,
            end,
            start_byte,
            end_byte,
        })
    }
}

/// Acceptable argument types for LAMMPS commands
#[derive(Debug)]
pub enum Argument {
    /// Expression from LAMMPS
    Int(isize),
    Float(f64),
    Bool(bool),
    ArgName(String), // TODO Rename to keyword arg?
    /// Variables within curly braces
    VarCurly(Ident),
    VarRound(expressions::Expression),
    String,
    Expression(expressions::Expression),
    // TODO Remove? Can't know if a group name until further on in the process???
    // Perhaps make it an identifier that then is decided to be either
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    /// TODO make this hold and `Ident` struct, from `crate::identifinder` module
    UnderscoreIdent(Ident),
}

impl Argument {
    fn from_node(
        node: &Node,
        // _cursor: &mut tree_sitter::TreeCursor,
        text: &[u8],
    ) -> Result<Self, String> {
        // TODO make these variants more complete.

        match node.child(0).unwrap().kind() {
            "int" => Ok(Self::Int(
                isize::from_str_radix(&node.child(0).unwrap().utf8_text(text).unwrap(), 10)
                    .unwrap(),
            )),
            "expression" => Ok(Self::Expression(
                // TODO get rid of this unwrap
                expressions::Expression::parse_expression(&node, text).unwrap(),
            )),
            "string" => Ok(Self::String),
            "group" => Ok(Self::Group),
            "underscore_ident" => Ok(Self::UnderscoreIdent(
                Ident::new(&node.child(0).unwrap(), text).map_err(|x| format!("{}", x))?,
            )),
            "var_curly" => Ok(Self::VarCurly(
                Ident::new(&node.child(0).unwrap().child(1).unwrap(), text)
                    .map_err(|x| format!("{}", x))?,
            )),
            "argname" => Ok(Self::ArgName(
                node.child(0).unwrap().utf8_text(text).unwrap().to_string(),
            )),
            x => Err(format!("Unknown argument type: {}", x)),
        }
    }
}

/// Commands that have a special form in the tree sitter grammar
/// TODO add arguments
/// TODO Add command location
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
            s => Err(format!("Unknon command: {s}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use tree_sitter::Parser;

    use super::ts_to_ast;

    fn setup_parser() -> Parser {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        parser
    }

    #[test]
    fn test_ast() {
        let mut parser = setup_parser();
        let source_bytes = include_bytes!("../../fix.lmp");
        let tree = parser.parse(source_bytes, None).unwrap();

        let ast = ts_to_ast(&tree, source_bytes);
        dbg!(ast);
        unimplemented!()
    }
}

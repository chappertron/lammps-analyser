//! Convert the treesiterr trees into an AST.
//! TODO Work out why the generic command nodes are causing errors:

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

    //}

    // todo!()
    //
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
    VarRound(Expression),
    String,
    Expression(Expression),
    // TODO Remove? Can't know if a group name until further on in the process???
    // Perhaps make it an identifier that then is decided to be either
    Group,
    /// A variable/fix/compute name prefixed by v_/f_/c_
    /// TODO make this hold and `Ident` struct, from `crate::identifinder` module
    UnderscoreIdent(Ident),
}

#[derive(Debug, Default, PartialEq, Clone)]
/// TODO convert into an enum
pub enum Expression {
    #[default]
    None,
    /// LAMMPS Identifier for a fix/compute/variable that is
    /// proceeded by f_/c_/v_ to indicate the type.
    UnderscoreIdent(Ident),
    /// An integer
    Int(isize),
    /// A float. Parsed as a 64-bit float
    Float(f64),
    /// A binary expression between two other expressions.
    BinaryOp(Box<Expression>, BinaryOp, Box<Expression>),
    // TODO MORE EXPRESSIONS
}

#[derive(Debug)]
enum ParseExprError {
    Utf8Error(std::str::Utf8Error),
    ParseIntError(std::num::ParseIntError),
    ParseFloatError(std::num::ParseFloatError),
}

impl From<std::num::ParseFloatError> for ParseExprError {
    fn from(v: std::num::ParseFloatError) -> Self {
        Self::ParseFloatError(v)
    }
}

impl From<std::num::ParseIntError> for ParseExprError {
    fn from(v: std::num::ParseIntError) -> Self {
        Self::ParseIntError(v)
    }
}

impl From<std::str::Utf8Error> for ParseExprError {
    fn from(v: std::str::Utf8Error) -> Self {
        Self::Utf8Error(v)
    }
}

impl Expression {
    /// TODO Handle Erros
    fn parse_expression(node: &Node<'_>, text: &[u8]) -> Result<Expression, ParseExprError> {
        dbg!(node);

        match node.kind() {
            "binary_op" => Ok(Self::BinaryOp(
                Box::new(Self::parse_expression(&node.child(0).unwrap(), text)?),
                // TODO Find a way to get the operator from teh TS symbol
                BinaryOp::try_from(node.child(1).unwrap().utf8_text(text)?).unwrap(),
                Box::new(Self::parse_expression(&node.child(2).unwrap(), text)?),
            )),
            "int" => Ok(Self::Int(isize::from_str_radix(
                &node.utf8_text(text)?,
                10,
            )?)),

            "float" => Ok(Self::Float(str::parse(&node.utf8_text(text)?)?)),

            // Just go into next level down
            "expression" => Ok(Self::parse_expression(&node.child(0).unwrap(), text)?),

            x => unimplemented!("Unknown expression type: {}", x),
        }
        // todo!()
    }
}
#[derive(Debug, Default, PartialEq, Eq, Clone, Copy)]
pub enum BinaryOp {
    #[default]
    Add,
    Subtract,
    Multiply,
    Divide,
    Power,
    Modulo,
    Equal,
    NotEqual,
    LessThan,
    LessThanOrEqual,
    GreaterThan,
    GreaterThanOrEqual,
    And,
    Or,
    Xor,
}
impl TryFrom<&str> for BinaryOp {
    type Error = String;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "+" => Ok(Self::Add),
            "-" => Ok(Self::Subtract),
            "*" => Ok(Self::Multiply),
            "/" => Ok(Self::Divide),
            "^" => Ok(Self::Power),
            "%" => Ok(Self::Modulo),
            "==" => Ok(Self::Equal),
            "!=" => Ok(Self::NotEqual),
            "<" => Ok(Self::LessThan),
            "<=" => Ok(Self::LessThanOrEqual),
            ">" => Ok(Self::GreaterThan),
            ">=" => Ok(Self::GreaterThanOrEqual),
            "&&" => Ok(Self::And),
            "||" => Ok(Self::Or),
            "^|" => Ok(Self::Xor),
            _ => Err(format!("Unknown binary operator: {}", value)),
        }
    }
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
                Expression::parse_expression(&node, text).unwrap(),
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

    use super::{ts_to_ast, BinaryOp, Expression};

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

    #[test]
    fn parse_expr() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 1+2";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();
        dbg!(command_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(command_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Int(1)),
                BinaryOp::Add,
                Box::new(Expression::Int(2))
            )
        );
    }

    #[test]
    fn parse_expr_floats() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 1.0+2.0";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();
        dbg!(command_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(command_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Float(1.0)),
                BinaryOp::Add,
                Box::new(Expression::Float(2.0))
            )
        );
    }

    #[test]
    fn parse_nested_expr() {
        let mut parser = setup_parser();
        let source_bytes = b"variable a equal 0.5*(1.0+2.0)";

        let tree = parser.parse(source_bytes, None).unwrap();

        // let ast = ts_to_ast(&tree, source_bytes);

        let command_node = tree.root_node().child(0).unwrap();
        dbg!(command_node.to_sexp());
        // Lots of tedium to parsing this...
        let expr_node = dbg!(command_node.child(0)).unwrap().child(3).unwrap();

        dbg!(expr_node.to_sexp());
        // assert_eq!(ast.commands.len(), 1);
        assert_eq!(
            // ast.commands[0],
            Expression::parse_expression(&expr_node, source_bytes.as_slice()).unwrap(),
            Expression::BinaryOp(
                Box::new(Expression::Float(0.5)),
                BinaryOp::Multiply,
                Box::new(Expression::BinaryOp(
                    Box::new(Expression::Float(1.0)),
                    BinaryOp::Add,
                    Box::new(Expression::Float(2.0))
                ))
            )
        );
    }
}

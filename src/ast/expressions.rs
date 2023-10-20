use tree_sitter::Node;

use crate::identifinder::Ident;

#[derive(Debug, Default, PartialEq, Clone)]
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
    /// An expression wrapped in brackets.
    /// TODO Perhaps just ignore brackets?
    Parens(Box<Expression>),
    // TODO MORE EXPRESSIONS
}

#[derive(Debug)]
pub(crate) enum ParseExprError {
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
    pub(crate) fn parse_expression(
        node: &Node<'_>,
        text: &[u8],
    ) -> Result<Expression, ParseExprError> {
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
            "parens" => Ok(Self::Parens(Box::new(Self::parse_expression(
                &node.child(1).unwrap(),
                text,
            )?))),

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

#[cfg(test)]
mod tests {
    use super::*;
    use tree_sitter::Parser;
    fn setup_parser() -> Parser {
        let mut parser = Parser::new();

        parser
            .set_language(tree_sitter_lammps::language())
            .expect("Could not load language");
        parser
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

        dbg!(tree.root_node().to_sexp());
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
                Box::new(Expression::Parens(Box::new(Expression::BinaryOp(
                    Box::new(Expression::Float(1.0)),
                    BinaryOp::Add,
                    Box::new(Expression::Float(2.0))
                ))))
            )
        );
    }
}

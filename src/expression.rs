use crate::gates::Parameter;

/// An enumeration describing arithmetic expressions that can be evaluated at
/// a later time.
pub enum Expression
{
    /// A direct value
    Value(Parameter),
    /// The sum of two expressions
    Sum(Box<Expression>, Box<Expression>),
    /// The difference between two expressions
    Difference(Box<Expression>, Box<Expression>),
    /// The product of two expressions
    Product(Box<Expression>, Box<Expression>),
    /// The quotient of two expressions
    Quotient(Box<Expression>, Box<Expression>),
    /// The negative of an expression
    Negative(Box<Expression>),
    /// One expression raised to the power of another
    Power(Box<Expression>, Box<Expression>),
    /// A function applied to an expression
    Function(String, Box<Expression>),
    /// A variable
    Variable(String)
}

impl Expression
{
    /// Create a new expression from fixed value `x`
    pub fn value<T>(x: T) -> Self
    where Parameter: From<T>
    {
        Expression::Value(Parameter::from(x))
    }

    /// Create a new expression, describing the sum of `left` and `right`
    pub fn sum(left: Expression, right: Expression) -> Self
    {
        Expression::Sum(Box::new(left), Box::new(right))
    }

    /// Create a new expression, describing the difference between `left` and `right`
    pub fn difference(left: Expression, right: Expression) -> Self
    {
        Expression::Difference(Box::new(left), Box::new(right))
    }

    /// Create a new expression, describing the product of `left` and `right`
    pub fn product(left: Expression, right: Expression) -> Self
    {
        Expression::Product(Box::new(left), Box::new(right))
    }

    /// Create a new expression, describing the `left` divided by `right`
    pub fn quotient(left: Expression, right: Expression) -> Self
    {
        Expression::Quotient(Box::new(left), Box::new(right))
    }

    /// Create a new expression, describing negative of `expr`
    pub fn negative(expr: Expression) -> Self
    {
        Expression::Negative(Box::new(expr))
    }

    /// Create a new expression, describing `left` raised to the power `right`
    pub fn power(left: Expression, right: Expression) -> Self
    {
        Expression::Power(Box::new(left), Box::new(right))
    }

    /// Create a new expression, describing the application of `function` on `arg`
    pub fn function(fname: String, arg: Expression) -> Self
    {
        Expression::Function(fname, Box::new(arg))
    }

    /// Create a new expression referencing an (as yet) unknown parameter `name`
    pub fn variable(name: String) -> Self
    {
        Expression::Variable(name)
    }

    /// Parse a real number
    ///
    /// Parse a real number. This can be either in the form of a real literal,
    /// an integer literal (which is converted to a real number), or the string
    /// "pi", which obviously indicates the number π. On success, the parsed
    /// number is returned, together with the remainder of the string to be
    /// parsed. On failure, a ParseError::InvalidArgument error is returned.
    fn parse_real_literal(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let real = regex::Regex::new(
            r"^\s*((?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?)"
        ).unwrap();
        let integer = regex::Regex::new(r"^\s*([1-9][0-9]*|0)").unwrap();
        let pi = regex::Regex::new(r"^\s*pi").unwrap();

        if let Some(captures) = real.captures(expr)
        {
            let m = captures.get(0).unwrap();
            if let Ok(nr) = captures[1].parse::<f64>()
            {
                Ok((Self::value(nr), &expr[m.end()..]))
            }
            else
            {
                // It's probably impossible to get here with the above regular
                // expression. Too large numbers are mapped to Inf.
                Err(crate::error::ParseError::InvalidArgument(String::from(expr)))
            }
        }
        else if let Some(captures) = integer.captures(expr)
        {
            let m = captures.get(0).unwrap();
            if let Ok(nr) = captures[1].parse::<u64>()
            {
                Ok((Self::value(nr as f64), &expr[m.end()..]))
            }
            else
            {
                Err(crate::error::ParseError::InvalidArgument(String::from(expr)))
            }
        }
        else if let Some(m) = pi.find(expr)
        {
            Ok((Self::value(::std::f64::consts::PI), &expr[m.end()..]))
        }
        else
        {
            Err(crate::error::ParseError::InvalidArgument(String::from(expr)))
        }
    }

    /// Parse a possibly parenthesized expression
    ///
    /// Parse an argument to a gate. This function expects either an expression
    /// in parentheses, or a literal real number. On success, the parsed
    /// number is returned, together with the remainder of the string to be
    /// parsed. On failure, a ParseError is returned.
    fn parse_parenthesized_expression(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let fun_open = regex::Regex::new(r"^\s*\(").unwrap();
        let fun_close = regex::Regex::new(r"^\s*\)").unwrap();
        if let Some(m) = fun_open.find(expr)
        {
            let (result, rest) = Self::parse_sum_expression(&expr[m.end()..])?;
            if let Some(m) = fun_close.find(rest)
            {
                Ok((result, &rest[m.end()..]))
            }
            else
            {
                Err(crate::error::ParseError::UnclosedParentheses(String::from(expr)))
            }
        }
        else
        {
            Self::parse_real_literal(expr)
        }
    }

    /// Parse a function call
    ///
    /// Parse an argument to a gate. If a function call is found, it is parsed
    /// in this function, otherwise control is passed on to
    /// `parse_parenthesized_expression()`. The functions that are recognised
    /// are `sin`, `cos`, `tan`, `exp`, `ln`, and `sqrt`. On success, the parsed
    /// number is returned, together with the remainder of the string to be
    /// parsed. On failure, a ParseError is returned.
    fn parse_function_expression(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let fun_open = regex::Regex::new(r"^\s*(sin|cos|tan|exp|ln|sqrt)\s*\(").unwrap();
        let fun_close = regex::Regex::new(r"^\s*\)").unwrap();
        if let Some(captures) = fun_open.captures(expr)
        {
            let m = captures.get(0).unwrap();
            let (arg, new_rest) = Self::parse_sum_expression(&expr[m.end()..])?;
            if let Some(m) = fun_close.find(new_rest)
            {
                Ok((Self::function(String::from(&captures[1]), arg), &new_rest[m.end()..]))
            }
            else
            {
                Err(crate::error::ParseError::UnclosedParentheses(String::from(expr)))
            }
        }
        else
        {
            Self::parse_parenthesized_expression(expr)
        }
    }

    /// Parse a power-raising expression
    ///
    /// Parse an argument to a gate, in the form of a number possibly raised
    /// to another number.  On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_power_expression(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let op = regex::Regex::new(r"^\s*\^").unwrap();
        let (left, rest) = Self::parse_function_expression(expr)?;
        if let Some(m) = op.find(rest)
        {
            let (right, new_rest) = Self::parse_power_expression(&rest[m.end()..])?;
            Ok((Self::power(left, right), new_rest))
        }
        else
        {
            Ok((left, rest))
        }
    }

    /// Parse a negated expression
    ///
    /// Parse an argument to a gate, in the form of a number that is zero or
    /// more times negated.  On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_negative_expression(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let op = regex::Regex::new(r"^\s*\-").unwrap();
        let mut rest = expr;
        let mut flip_sign = false;
        while let Some(m) = op.find(rest)
        {
            rest = &rest[m.end()..];
            flip_sign = !flip_sign;
        }

        let (result, new_rest) = Self::parse_power_expression(rest)?;
        if flip_sign
        {
            Ok((Self::negative(result), new_rest))
        }
        else
        {
            Ok((result, new_rest))
        }
    }

    /// Parse a product expression
    ///
    /// Parse an argument to a gate, in the form of a product or quotient of one
    /// or more expressions. On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_product_expression(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let op = regex::Regex::new(r"^\s*([*/])").unwrap();
        let (mut left, mut rest) = Self::parse_negative_expression(expr)?;
        while let Some(captures) = op.captures(rest)
        {
            let m = captures.get(0).unwrap();
            let (right, new_rest) = Self::parse_negative_expression(&rest[m.end()..])?;
            if &captures[1] == "*"
            {
                left = Self::product(left, right);
            }
            else
            {
                left = Self::quotient(left, right);
            }
            rest = new_rest;
        }

        Ok((left, rest))
    }

    /// Parse a sum expression
    ///
    /// Parse an expression, in the form of a sum or difference of one
    /// or more expressions. On success, the parsed number is returned, together
    /// with the remainder of the string to be parsed. On failure, a ParseError
    /// is returned.
    fn parse_sum_expression(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        let op = regex::Regex::new(r"^\s*([-+])").unwrap();
        let (mut left, mut rest) = Self::parse_product_expression(expr)?;
        while let Some(captures) = op.captures(rest)
        {
            let m = captures.get(0).unwrap();
            let (right, new_rest) = Self::parse_product_expression(&rest[m.end()..])?;
            if &captures[1] == "+"
            {
                left = Self::sum(left, right);
            }
            else
            {
                left = Self::difference(left, right);
            }
            rest = new_rest;
        }

        Ok((left, rest))
    }

    /// Parse an expression
    ///
    /// Parse an arithmetic expression. On success, the parsed expression is
    /// returned, together with the remainder of the string to be parsed. On
    /// failure, a ParseError is returned.
    pub fn parse(expr: &str) -> crate::error::ParseResult<(Self, &str)>
    {
        Self::parse_sum_expression(expr)
    }

    /// Evaluate an expression
    ///
    /// Evaluate this expression and return the result. In case the expression
    /// uses any variables, return an UnknownVariable error.
    pub fn eval(&self) -> crate::error::Result<f64>
    {
        self.eval_with_parameters(&vec![])
    }

    /// Evaluate an expression
    ///
    /// Evaluate this expression, substituting the values in `param_map` for the
    /// variables in the expression, and return the result. In case a variable
    /// in the expression is not in `param_map`, an UnknownVariable error is
    /// returned. In case a function name outside the known set (`sin`, `cos`,
    /// '`tan`, `exp`, `ln`, and `sqrt`) is used, an UnknownFunction error
    /// is returned.
    pub fn eval_with_parameters(&self, param_map: &Vec<(&String, f64)>)
        -> crate::error::Result<f64>
    {
        match *self
        {
            Expression::Value(ref x) => Ok(x.value()),
            Expression::Sum(ref left, ref right) => {
                Ok(left.eval_with_parameters(param_map)?
                    + right.eval_with_parameters(&param_map)?)
            },
            Expression::Difference(ref left, ref right) => {
                Ok(left.eval_with_parameters(param_map)?
                    - right.eval_with_parameters(&param_map)?)
            },
            Expression::Product(ref left, ref right) => {
                Ok(left.eval_with_parameters(param_map)?
                    * right.eval_with_parameters(&param_map)?)
            },
            Expression::Quotient(ref left, ref right) => {
                Ok(left.eval_with_parameters(param_map)?
                    / right.eval_with_parameters(&param_map)?)
            },
            Expression::Negative(ref expr) => {
                Ok(-expr.eval_with_parameters(param_map)?)
            },
            Expression::Power(ref left, ref right) => {
                Ok(left.eval_with_parameters(param_map)?
                    .powf(right.eval_with_parameters(&param_map)?))
            },
            Expression::Function(ref fname, ref arg) => {
                let x = arg.eval_with_parameters(param_map)?;
                match &**fname
                {
                    "sin"  => Ok(x.sin()),
                    "cos"  => Ok(x.cos()),
                    "tan"  => Ok(x.tan()),
                    "exp"  => Ok(x.exp()),
                    "ln"   => Ok(x.ln()),
                    "sqrt" => Ok(x.sqrt()),
                    _      => {
                        Err(crate::error::Error::UnknownFunction(fname.clone()))
                    }

                }
            },
            Expression::Variable(ref name) => {
                for (pname, pvalue) in param_map
                {
                    if pname == &name
                    {
                        return Ok(*pvalue);
                    }
                }
                Err(crate::error::Error::UnknownVariable(name.clone()))
            }
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::Expression;
    use crate::gates::Parameter;
    use num_traits::Float;

    fn is_near(a: f64, b: f64) -> bool
    {
        let scale = a.abs().max(b.abs());
        (a-b).abs() <= scale * f64::epsilon()
    }

    #[test]
    fn test_constructors()
    {
        assert!(matches!(Expression::value(1.12), Expression::Value(Parameter::Direct(_))));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        assert!(matches!(Expression::sum(a, b), Expression::Sum(_, _)));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        assert!(matches!(Expression::difference(a, b), Expression::Difference(_, _)));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        assert!(matches!(Expression::product(a, b), Expression::Product(_, _)));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        assert!(matches!(Expression::quotient(a, b), Expression::Quotient(_, _)));

        let a = Expression::value(1.12);
        assert!(matches!(Expression::negative(a), Expression::Negative(_)));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        assert!(matches!(Expression::power(a, b), Expression::Power(_, _)));

        let a = Expression::value(1.12);
        assert!(matches!(
            Expression::function(String::from("sin"), a),
            Expression::Function(_, _)
        ));

        assert!(matches!(
            Expression::variable(String::from("x")),
            Expression::Variable(_)
        ));
    }

    #[test]
    fn test_eval()
    {
        let a = Expression::value(-0.23);
        assert!(is_near(a.eval().unwrap(), -0.23));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        let c = Expression::sum(a, b);
        assert!(is_near(c.eval().unwrap(), 5.08));

        let a = Expression::value(1.12);
        let b = Expression::value(3.96);
        let c = Expression::difference(a, b);
        assert!(is_near(c.eval().unwrap(), -2.84));

        let a = Expression::value(0.11);
        let b = Expression::value(3.14);
        let c = Expression::product(a, b);
        assert!(is_near(c.eval().unwrap(), 0.3454));

        let a = Expression::value(3.14);
        let b = Expression::value(-0.5);
        let c = Expression::quotient(a, b);
        assert!(is_near(c.eval().unwrap(), -6.28));

        let a = Expression::value(3.14);
        let b = Expression::negative(a);
        assert!(is_near(b.eval().unwrap(), -3.14));

        let a = Expression::value(3.14);
        let b = Expression::value(0.4);
        let c = Expression::power(a, b);
        assert!(is_near(c.eval().unwrap(), 1.580417606406723));

        let a = Expression::variable(String::from("a"));
        assert!(matches!(a.eval(), Err(crate::error::Error::UnknownVariable(_))));
    }

    #[test]
    fn test_functions()
    {
        let x = Expression::value(1.23);
        let a = Expression::function(String::from("sin"), x);
        assert!(is_near(a.eval().unwrap(), 0.9424888019316975));

        let x = Expression::value(1.23);
        let a = Expression::function(String::from("cos"), x);
        assert!(is_near(a.eval().unwrap(), 0.3342377271245026));

        let x = Expression::value(1.23);
        let a = Expression::function(String::from("tan"), x);
        assert!(is_near(a.eval().unwrap(), 2.819815734268152));

        let x = Expression::value(1.23);
        let a = Expression::function(String::from("exp"), x);
        assert!(is_near(a.eval().unwrap(), 3.4212295362896734));

        let x = Expression::value(1.23);
        let a = Expression::function(String::from("ln"), x);
        assert!(is_near(a.eval().unwrap(), 0.20701416938432612));

        let x = Expression::value(1.23);
        let a = Expression::function(String::from("sqrt"), x);
        assert!(is_near(a.eval().unwrap(), 1.1090536506409416));
    }

    #[test]
    fn test_parse()
    {
        let res = Expression::parse("1 + 2 * 3");
        assert!(matches!(res, Ok(_)));
        let (x, rest) = res.unwrap();
        assert!(is_near(x.eval().unwrap(), 7.0));
        assert_eq!(rest, "");

        let res = Expression::parse("1/2 - (1+4)");
        assert!(matches!(res, Ok(_)));
        let (x, rest) = res.unwrap();
        assert!(is_near(x.eval().unwrap(), -4.5));
        assert_eq!(rest, "");

        let res = Expression::parse("1/2 - (1+4)");
        assert!(matches!(res, Ok(_)));
        let (x, rest) = res.unwrap();
        assert!(is_near(x.eval().unwrap(), -4.5));
        assert_eq!(rest, "");

        let res = Expression::parse("sin(1/2)");
        assert!(matches!(res, Ok(_)));
        let (x, rest) = res.unwrap();
        assert!(is_near(x.eval().unwrap(), 0.479425538604203));
        assert_eq!(rest, "");
    }
}

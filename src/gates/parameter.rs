// Copyright 2019 Q1t BV
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/// Enumeration for parameters (usually angles) to a gate. This can either be a
/// direct floating points value, or a reference to one. In the latter case,
/// the value of the parameter can be changed outside the circuit, which can
/// be useful to e.g. to reexecute a circuit with the end state of the last
/// excution, but with a different value for a parameter. For reference
/// parameters, a name is associated with the parameter for the description,
/// since its value can change.
pub enum Parameter
{
    /// Direct value
    Direct(f64),
    /// Reference value, mutable outside the circuit, with its name
    Reference(::std::rc::Rc<::std::cell::RefCell<f64>>, String),
    /// Reference parameter from external code
    FFIRef(*const f64)
}

impl Parameter
{
    /// Create a new reference parameter from `cell` with name  `name`.
    pub fn from_refcell(cell: &::std::rc::Rc<::std::cell::RefCell<f64>>, name: &str) -> Self
    {
        Parameter::Reference(cell.clone(), String::from(name))
    }

    /// Return the current value of the parameter
    pub fn value(&self) -> f64
    {
        match *self
        {
            Parameter::Direct(p) => p,
            Parameter::Reference(ref p, _) => *p.borrow(),
            Parameter::FFIRef(p) => unsafe { *p }
        }
    }
}

impl From<f64> for Parameter
{
    fn from(val: f64) -> Self
    {
        Parameter::Direct(val)
    }
}

impl ::std::fmt::Display for Parameter
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> std::fmt::Result
    {
        match *self
        {
            Parameter::Direct(p) => p.fmt(f),
            Parameter::Reference(_, ref name) => write!(f, "{}", name),
            Parameter::FFIRef(ptr) => {
                let p = unsafe { *ptr };
                p.fmt(f)
            }
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::Parameter;

    #[test]
    fn test_from()
    {
        let p = Parameter::from(1.23);
        assert!(matches!(p, Parameter::Direct(_)));
        assert_eq!(p.value(), 1.23);

        let p = Parameter::from(-9.23);
        assert!(matches!(p, Parameter::Direct(_)));
        assert_eq!(p.value(), -9.23);
    }

    #[test]
    fn test_value()
    {
        let cell = ::std::rc::Rc::new(::std::cell::RefCell::new(0.12));
        let x = 3.14_f64;

        let p0 = Parameter::from(x);
        let p1 = Parameter::from_refcell(&cell, "x");
        assert_eq!(p0.value(), 3.14);
        assert_eq!(p1.value(), 0.12);

        *cell.borrow_mut() = 6.24;
        assert_eq!(p0.value(), 3.14);
        assert_eq!(p1.value(), 6.24);
    }

    #[test]
    fn test_display()
    {
        let cell = ::std::rc::Rc::new(::std::cell::RefCell::new(0.12));
        let x = 3.14_f64;

        let p0 = Parameter::from(x);
        let p1 = Parameter::from_refcell(&cell, "x");
        let p2 = Parameter::from_refcell(&cell, "longname");

        assert_eq!(format!("{}", p0), String::from("3.14"));
        assert_eq!(format!("{}", p1), String::from("x"));
        assert_eq!(format!("{}", p2), String::from("longname"));
        assert_eq!(format!("{:.4}", p0), String::from("3.1400"));
        assert_eq!(format!("{:.4}", p1), String::from("x"));
        assert_eq!(format!("{:.4}", p2), String::from("longname"));
    }
}

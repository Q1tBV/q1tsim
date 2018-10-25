extern crate num_complex;

use cmatrix;
use gates;

/// The `V` gate
///
/// The `V` gate is the square root of the `X` gate.
pub struct V
{
}

impl V
{
    /// Create a new `V` gate.
    pub fn new() -> Self
    {
        V { }
    }
}

impl gates::Gate for V
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "V"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let h = 0.5 * cmatrix::COMPLEX_ONE;
        let hi = 0.5 * cmatrix::COMPLEX_I;
        array![[h+hi, h-hi], [h-hi, h+hi]]
    }

    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("u3(pi/2, -pi/2, pi/2) {}", bit_names[bits[0]])
    }

    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        let name = &bit_names[bits[0]];
        format!("sdag {}\nh {}\nsdag {}", name, name, name)
    }
}

/// Conjugate of `V` gate.
///
/// The `V`<sup>`†`</sup> is the conjugate of the `V` gate.
pub struct Vdg
{
}

impl Vdg
{
    /// Create a new `V`<sup>`†`</sup> gate.
    pub fn new() -> Self
    {
        Vdg { }
    }
}

impl gates::Gate for Vdg
{
    fn cost(&self) -> f64
    {
        gates::U3::cost()
    }

    fn description(&self) -> &str
    {
        "V†"
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let h = 0.5 * cmatrix::COMPLEX_ONE;
        let hi = 0.5 * cmatrix::COMPLEX_I;
        array![[h-hi, h+hi], [h+hi, h-hi]]
    }

    fn open_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        format!("u3(pi/2, pi/2, -pi/2) {}", bit_names[bits[0]])
    }

    fn c_qasm(&self, bit_names: &[String], bits: &[usize]) -> String
    {
        let name = &bit_names[bits[0]];
        format!("s {}\nh {}\ns {}", name, name, name)
    }
}

#[cfg(test)]
mod tests
{
    use gates::{gate_test, Gate, V, Vdg};
    use cmatrix;

    #[test]
    fn test_description()
    {
        let gate = V::new();
        assert_eq!(gate.description(), "V");
        let gate = Vdg::new();
        assert_eq!(gate.description(), "V†");
    }

    #[test]
    fn test_matrix()
    {
        let h = 0.5*cmatrix::COMPLEX_ONE;
        let hi = 0.5*cmatrix::COMPLEX_I;

        let gate = V::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[h+hi, h-hi], [h-hi, h+hi]]);

        let gate = Vdg::new();
        assert_complex_matrix_eq!(gate.matrix(), array![[h-hi, h+hi], [h+hi, h-hi]]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let i = cmatrix::COMPLEX_I;
        let x = cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * o;
        let hi = 0.5 * i;

        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [h+hi, h-hi, x,  x*i],
            [h-hi, h+hi, x, -x*i]
        ];
        gate_test(V::new(), &mut state, &result);

        let mut state = array![
            [o, z, x,  x],
            [z, o, x, -x]
        ];
        let result = array![
            [h-hi, h+hi, x, -x*i],
            [h+hi, h-hi, x,  x*i]
        ];
        gate_test(Vdg::new(), &mut state, &result);
    }

    #[test]
    fn test_open_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = V::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "u3(pi/2, -pi/2, pi/2) qb");
        let qasm = Vdg::new().open_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "u3(pi/2, pi/2, -pi/2) qb");
    }

    #[test]
    fn test_c_qasm()
    {
        let bit_names = [String::from("qb")];
        let qasm = V::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "sdag qb\nh qb\nsdag qb");
        let qasm = Vdg::new().c_qasm(&bit_names, &[0]);
        assert_eq!(qasm, "s qb\nh qb\ns qb");
    }
}

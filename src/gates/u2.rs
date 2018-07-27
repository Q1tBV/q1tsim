extern crate num_complex;

use cmatrix;
use gates;

/// U<sub>2</sub> gate.
///
/// The `U`<sub>`2`</sub>`(ϕ, λ)` gate transforms a qubit by the matrix
/// ```text
///           ┌                      ┐
///           │      1       -exp(iλ)│
/// 1/sqrt(2) │                      │
///           │exp(iϕ)    exp(i(λ+ϕ))│
///           └                      ┘
/// ```
pub struct U2
{
    phi: f64,
    lambda: f64,
    desc: String
}

impl U2
{
    /// Create a new `U`<sub>`2`</sub> gate.
    pub fn new(phi: f64, lambda: f64) -> Self
    {
        U2 { phi: phi, lambda: lambda, desc: format!("U2({:.4}, {:.4})", phi, lambda) }
    }
}

impl gates::Gate for U2
{
    fn description(&self) -> &str
    {
        &self.desc
    }

    fn nr_affected_bits(&self) -> usize
    {
        1
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        array![[ num_complex::Complex::new(x, 0.0),
                -num_complex::Complex::from_polar(&x, &self.lambda)],
               [ num_complex::Complex::from_polar(&x, &self.phi),
                 num_complex::Complex::from_polar(&x, &(self.phi+self.lambda))]]
    }
}

#[cfg(test)]
mod tests
{
    extern crate num_complex;

    use gates::{gate_test, Gate, U2};
    use cmatrix;
    use self::num_complex::Complex;

    #[test]
    fn test_description()
    {
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_eq!(gate.description(), "U2(0.7854, 0.6931)");
    }

    #[test]
    fn test_matrix()
    {
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_complex_matrix_eq!(gate.matrix(), array![
            [Complex::new(0.7071067811865476, 0.0), Complex::new(-0.5439340435069544, -0.4518138513969824)],
            [Complex::new(               0.5, 0.5), Complex::new(0.06513881252516862,  0.7041000888388035)]
        ]);
    }

    #[test]
    fn test_apply()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = array![[x, x], [z, -x], [x, z], [z, z]];
        let result = array![
            [    o-i,        o],
            [      z,       -o],
            [2.0*i*x,  (o+i)*x],
            [      z, -(o+i)*x]
        ] * (0.5 * o);
        let gate = U2::new(::std::f64::consts::FRAC_PI_4, ::std::f64::consts::FRAC_PI_2);
        gate_test(gate, &mut state, &result);
    }
}

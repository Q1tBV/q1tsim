extern crate num_complex;
extern crate rulinalg;

use cmatrix;
use gates;

/// U<sub>3</sub> gate.
///
/// The `U`<sub>`3`</sub>`(θ, ϕ, λ)` gate is the univeral single-qubit
/// transformation, all unary gates can be written in terms of `U`<sub>`3`</sub>.
/// It transforms a qubit by the matrix
/// ```text
/// ┌                                     ┐
/// │       cos(θ/2)      -exp(iλ)sin(θ/2)│
/// │                                     │
/// │exp(iϕ)sin(θ/2)   exp(i(λ+ϕ))cos(θ/2)│
/// └                                     ┘
/// ```
pub struct U3
{
    theta: f64,
    phi: f64,
    lambda: f64,
    desc: String
}

impl U3
{
    /// Create a new `U`<sub>`3`</sub> gate.
    pub fn new(theta: f64, phi: f64, lambda: f64) -> Self
    {
        U3
        {
            theta: theta,
            phi: phi,
            lambda: lambda,
            desc: format!("U3({:.4}, {:.4}, {:.4})", theta, phi, lambda)
        }
    }
}

impl gates::Gate for U3
{
    fn description(&self) -> &str
    {
        &self.desc
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        let (c, s) = (self.theta.cos(), self.theta.sin());
        cmatrix::CMatrix::new(2, 2, vec![
            num_complex::Complex::new(c, 0.0),
            -num_complex::Complex::from_polar(&s, &self.lambda),
            num_complex::Complex::from_polar(&s, &self.phi),
            num_complex::Complex::from_polar(&c, &(self.phi+self.lambda))
        ])
    }
}

impl gates::UnaryGate for U3 { /* use default implementation */ }

#[cfg(test)]
mod tests
{
    extern crate num_complex;

    use gates::{Gate, UnaryGate};
    use gates::U3;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;
    use self::num_complex::Complex;

    #[test]
    fn test_description()
    {
        let gate = U3::new(0.17, ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_eq!(gate.description(), "U3(0.1700, 0.7854, 0.6931)");
    }

    #[test]
    fn test_matrix()
    {
        let gate = U3::new(0.32, ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::LN_2);
        assert_complex_matrix_eq!(gate.matrix().as_ref(), matrix![
            Complex::new(0.9492354180824408,                0.0), Complex::new(-0.2419768354941858, -0.2009958510568650);
            Complex::new(0.2224321481461860, 0.2224321481461860), Complex::new(0.08744374907134395,  0.9451991693238490)
        ]);
    }

    #[test]
    fn test_apply_unary()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let o = cmatrix::COMPLEX_ONE;
        let x = cmatrix::COMPLEX_HSQRT2;
        let i = cmatrix::COMPLEX_I;
        let mut state = cmatrix::CMatrix::new(4, 2, vec![x, x, z, -x, x, z, z, z]);

        let gate = U3::new(3.0*::std::f64::consts::FRAC_PI_4,
            ::std::f64::consts::FRAC_PI_4, ::std::f64::consts::FRAC_PI_2);

        gate.apply_unary(state.as_mut());
        assert_complex_matrix_eq!(state.as_ref(), matrix![
               -o-i,       -o;
                  z,        o;
              2.0*x,  (o+i)*x;
                  z, -(o+i)*x
        ] * (0.5 * o));
    }
}

extern crate rand;

use cmatrix;

/// Single quantum bit.
///
/// Struct Qubit represents a single quantum bit. It is a superposition of
/// the two basis states, &alpha;|0&rang; + &beta;|1&rang;, where &alpha; and
/// &beta; are complex numbers satisfying
/// |&alpha;|<sup>2</sup> + |&beta;|<sup>2</sup> = 1.
pub struct Qubit
{
    c: cmatrix::CMatrix
}

impl Qubit
{
    /// Create a new qubit, initialized to |0&rang;.
    pub fn new() -> Self
    {
        Qubit { c: cmatrix::CMatrix::new_real(2, 1, vec![1.0, 0.0]) }
    }

    /// Create a new qubit from coefficients.
    ///
    /// Create a new qubit initialized to (`ra`+`ia`i)|0&rang; + (`rb`+`ib`i)|1&rang;.
    /// The coefficients are scaled such that the qubit is normalized.
    pub fn from_coefs(ra: f64, ia: f64, rb: f64, ib: f64) -> Self
    {
        let c = if ia == 0.0 && ib == 0.0
            {
                let norm = ra.hypot(rb);
                cmatrix::CMatrix::new_real(2, 1, vec![ra/norm, rb/norm])
            }
            else if ra == 0.0 && rb == 0.0
            {
                let norm = ia.hypot(ib);
                cmatrix::CMatrix::new_imag(2, 1, vec![ia/norm, ib/norm])
            }
            else
            {
                let norm = (ra*ra + ia*ia + rb*rb + ib*ib).sqrt();
                cmatrix::CMatrix::new_complex(2, 1, vec![ra/norm, rb/norm], vec![ia/norm, ib/norm])
            };
        Qubit { c: c }
    }

    /// Measure the qubit.
    ///
    /// Perform a measurement on the qubit. If the qubit is in the state
    /// &alpha;|0&rang; + &beta;|1&rang;, a value of 0 is returned with
    /// probability |&alpha;|<sup>2</sup>, and a value of 1 with probability
    /// |&beta;|<sup>2</sup> = 1 - |&alpha;|<sup>2</sup>.
    pub fn measure(&self) -> u8
    {
        let weights = self.c.abs_sq();
        let x: f64 = rand::random();
        (x >= weights[[0,0]]) as u8
    }
}

#[cfg(test)]
mod tests
{
    use qubit::Qubit;

    #[test]
    fn test_measure()
    {
        let b = Qubit::new();
        assert_eq!(b.measure(), 0);

        let b = Qubit::from_coefs(0.0, 0.0, 0.0, 1.0);
        assert_eq!(b.measure(), 1);
    }
}

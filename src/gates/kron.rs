use cmatrix;
use gates;

/// Gate describing the Kronecker product of two other gates operating on
/// different qubits.
pub struct Kron<G0, G1>
{
    g0: G0,
    g1: G1,
    desc: String
}

impl<G0, G1> Kron<G0, G1>
where G0: gates::Gate, G1: gates::Gate
{
    /// Create a new Kronecker product gate `g1` ⊗ `g2`.
    pub fn new(g0: G0, g1: G1) -> Self
    {
        let desc = format!("{}⊗{}", g0.description(), g1.description());
        Kron { g0: g0, g1: g1, desc: desc }
    }
}

impl<G0, G1> gates::Gate for Kron<G0, G1>
where G0: gates::Gate, G1: gates::Gate
{
    fn description(&self) -> &str
    {
        &self.desc
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        self.g0.matrix().kron(&self.g1.matrix())
    }
}

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::Hadamard;
    use gates::Identity;
    use gates::Kron;
    use cmatrix;
    use rulinalg::matrix::BaseMatrix;

    #[test]
    fn test_description()
    {
        let ih = Kron::new(Identity::new(), Hadamard::new());
        assert_eq!(ih.description(), "I⊗H");
        let hh = Kron::new(Hadamard::new(), Hadamard::new());
        assert_eq!(hh.description(), "H⊗H");
        let hih = Kron::new(Hadamard::new(), Kron::new(Identity::new(), Hadamard::new()));
        assert_eq!(hih.description(), "H⊗I⊗H");
    }

    #[test]
    fn test_matrix()
    {
        let z = cmatrix::COMPLEX_ZERO;
        let s = cmatrix::COMPLEX_HSQRT2;
        let h = 0.5 * cmatrix::COMPLEX_ONE;

        let ih = Kron::new(Identity::new(), Hadamard::new());
        assert_complex_matrix_eq!(ih.matrix(), matrix![
            s,  s, z,  z;
            s, -s, z,  z;
            z,  z, s,  s;
            z,  z, s, -s
        ]);

        let hh = Kron::new(Hadamard::new(), Hadamard::new());
        assert_complex_matrix_eq!(hh.matrix(), matrix![
            h,  h,  h,  h;
            h, -h,  h, -h;
            h,  h, -h, -h;
            h, -h, -h,  h
        ]);

        let hih = Kron::new(Hadamard::new(), Kron::new(Identity::new(), Hadamard::new()));
        assert_complex_matrix_eq!(hih.matrix(), matrix![
            h,  h,  z,  z,  h,  h,  z,  z;
            h, -h,  z,  z,  h, -h,  z,  z;
            z,  z,  h,  h,  z,  z,  h,  h;
            z,  z,  h, -h,  z,  z,  h, -h;
            h,  h,  z,  z, -h, -h,  z,  z;
            h, -h,  z,  z, -h,  h,  z,  z;
            z,  z,  h,  h,  z,  z, -h, -h;
            z,  z,  h, -h,  z,  z, -h,  h
        ]);
    }
}

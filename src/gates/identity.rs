use cmatrix;
use gates;

/// The identity gate
///
/// The identity gate leaves the qubits on which it acts unchanged.
pub struct Identity
{
}

impl Identity
{
    /// Create a new identity gate.
    pub fn new() -> Self
    {
        Identity { }
    }
}

impl gates::Gate for Identity
{
    fn description(&self) -> &str
    {
        "I"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::CMatrix::eye(2, 2)
    }
}

#[cfg(test)]
mod tests
{
    use gates::Gate;
    use gates::Identity;

    #[test]
    fn test_description()
    {
        let i = Identity::new();
        assert_eq!(i.description(), "I");
    }

    #[test]
    fn test_matrix()
    {
        let i = Identity::new();
        assert_eq!(i.matrix().real().data(), &vec![1.0, 0.0, 0.0, 1.0]);
        assert_eq!(i.matrix().imag().data(), &vec![0.0; 4]);
    }
}

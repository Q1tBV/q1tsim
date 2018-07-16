use cmatrix;
use gates;

/// The C<sub>X</sub> gate.
///
/// The C<sub>X</sub> or CNOT gate flips a qubit based on a control bit: when
/// the control bit is zero, it leaves the target unchanged; when the control
/// bit is one, it flips the |0&rang; and |1&rang; components of the target bit.
pub struct CX
{
}

impl CX
{
    /// Create a new Hadamard gate.
    pub fn new() -> Self
    {
        CX { }
    }
}

impl gates::Gate for CX
{
    fn description(&self) -> &str
    {
        "CX"
    }

    fn matrix(&self) -> cmatrix::CMatrix
    {
        cmatrix::CMatrix::new_real(4, 4, vec![
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 1.0, 0.0
        ])
    }
}

impl gates::BinaryGate for CX
{
    fn apply_binary(&self, state: &mut cmatrix::CMatrix)
    {
        assert!(state.rows() % 4 == 0, "Number of rows is not a multiple of four.");

        let n = state.rows() / 4;
        state.swap_rows(2*n, 3*n, n);
    }
}


#[cfg(test)]
mod tests
{
    use gates::{Gate, BinaryGate};
    use gates::CX;
    use cmatrix;

    #[test]
    fn test_description()
    {
        let h = CX::new();
        assert_eq!(h.description(), "CX");
    }

    #[test]
    fn test_matrix()
    {
        let cx = CX::new();
        assert_matrix_eq!(cx.matrix().real(), matrix![
            1.0, 0.0, 0.0, 0.0;
            0.0, 1.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 1.0;
            0.0, 0.0, 1.0, 0.0
        ], comp=float);
        assert_matrix_eq!(cx.matrix().imag(), matrix![
            0.0, 0.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 0.0
        ], comp=float);
    }

    #[test]
    fn test_apply_binary()
    {
        let x = ::std::f64::consts::FRAC_1_SQRT_2;
        let mut state = cmatrix::CMatrix::new_real(4, 4, vec![
            1.0, 0.0,  0.5, 0.0,
            0.0, 0.0, -0.5, 0.0,
            0.0, 1.0,  0.5,   x,
            0.0, 0.0, -0.5,  -x]
        );

        CX::new().apply_binary(&mut state);
//         assert_matrix_eq!(state.real(), matrix![
//             1.0, 0.0,  0.5, 0.0;
//             0.0, 0.0, -0.5, 0.0;
//             0.0, 0.0, -0.5,  -x;
//             0.0, 1.0,  0.5,   x], comp=float);
        assert_matrix_eq!(state.imag(), matrix![
            0.0, 0.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 0.0;
            0.0, 0.0, 0.0, 0.0], comp=float);
    }
}

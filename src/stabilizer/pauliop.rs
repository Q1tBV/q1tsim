/// Enumeration for the Pauli spin operations
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum PauliOp
{
    /// Identity
    I,
    /// Pauli Z
    Z,
    /// Pauli X
    X,
    /// Pauli Y
    Y
}

impl PauliOp
{
    /// Create a new PauliOp from its bit representation
    pub fn from_bits(bits: u64) -> Self
    {
        match bits & 0x03
        {
            0 => PauliOp::I,
            1 => PauliOp::Z,
            2 => PauliOp::X,
            3 => PauliOp::Y,
            _ => unreachable!()
        }
    }

    /// Convert a PauliOp to its bit representation
    pub fn to_bits(self) -> u64
    {
        match self
        {
            PauliOp::I => 0,
            PauliOp::Z => 1,
            PauliOp::X => 2,
            PauliOp::Y => 3
        }
    }
}

impl ::std::fmt::Display for PauliOp
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result
    {
        match *self
        {
            PauliOp::I => write!(f, "I"),
            PauliOp::Z => write!(f, "Z"),
            PauliOp::X => write!(f, "X"),
            PauliOp::Y => write!(f, "Y"),
        }
    }
}

#[cfg(test)]
mod tests
{
    use super::PauliOp;

    #[test]
    fn test_from_bits()
    {
        assert_eq!(PauliOp::from_bits(0), PauliOp::I);
        assert_eq!(PauliOp::from_bits(1), PauliOp::Z);
        assert_eq!(PauliOp::from_bits(2), PauliOp::X);
        assert_eq!(PauliOp::from_bits(3), PauliOp::Y);
        assert_eq!(PauliOp::from_bits(4), PauliOp::I);
        assert_eq!(PauliOp::from_bits(0xff), PauliOp::Y);
    }

    #[test]
    fn test_to_bits()
    {
        assert_eq!(PauliOp::to_bits(PauliOp::I), 0);
        assert_eq!(PauliOp::to_bits(PauliOp::Z), 1);
        assert_eq!(PauliOp::to_bits(PauliOp::X), 2);
        assert_eq!(PauliOp::to_bits(PauliOp::Y), 3);
    }

    #[test]
    fn test_display()
    {
        assert_eq!(format!("{}", PauliOp::I), "I");
        assert_eq!(format!("{}", PauliOp::Z), "Z");
        assert_eq!(format!("{}", PauliOp::X), "X");
        assert_eq!(format!("{}", PauliOp::Y), "Y");
    }
}

#[macro_use] extern crate rulinalg;

pub mod cmatrix;
pub mod gates;

#[cfg(test)]
mod tests
{
    use gates;
    use gates::Gate;

    #[test]
    fn test_description()
    {
        let i = gates::Identity::new();
        assert_eq!(i.description(), "I");
        let h = gates::Hadamard::new();
        assert_eq!(h.description(), "H");
    }

    #[test]
    fn test_matrix()
    {
        let i = gates::Identity::new();
        assert_eq!(i.matrix().real(), Some(&matrix![1.0, 0.0; 0.0, 1.0]));
        assert_eq!(i.matrix().imag(), None);
        assert_eq!(i.matrix().real_expand().data(), &vec![1.0, 0.0, 0.0, 1.0]);
        assert_eq!(i.matrix().imag_expand().data(), &vec![0.0; 4]);

        let h = gates::Hadamard::new();
        let s = ::std::f64::consts::FRAC_1_SQRT_2;
        assert_eq!(h.matrix().real(), Some(&matrix![s, s; s, -s]));
        assert_eq!(h.matrix().imag(), None);
        assert_eq!(h.matrix().real_expand().data(), &vec![s, s, s, -s]);
        assert_eq!(h.matrix().imag_expand().data(), &vec![0.0; 4]);
    }
}

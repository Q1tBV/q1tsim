pub mod gates;

#[cfg(test)]
mod tests {
    use gates;
    use gates::Gate;

    #[test]
    fn test_description()
    {
        let h = gates::Hadamard::new();
        assert_eq!(h.description(), "H");
    }
}

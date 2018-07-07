mod hadamard;

pub trait Gate
{
    fn description(&self) -> &str;
}

pub use gates::hadamard::Hadamard;

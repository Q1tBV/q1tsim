use gates;

pub struct Hadamard
{
}

impl Hadamard
{
    pub fn new() -> Self
    {
        Hadamard { }
    }
}

impl gates::Gate for Hadamard
{
    fn description(&self) -> &str
    {
        return "H";
    }
}

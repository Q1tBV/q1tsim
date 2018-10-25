use gates;

/// Custom gate definition.
///
/// Struct CustomDef provides is used to create user-defined gates that are
/// made out of a sequence of more primitive gates.
pub struct CustomDef
{
    // The name of the gate
    name: String,
    // The number of gates on which this gate operates
    nr_bits: usize,
    // The number of parameters to this gate
    nr_params: usize,
    // The operations making up the gate
    ops: Vec<CustomDefOp>
}

struct GateBuilder<G> where G: ?Sized
{
    param_idxs: Vec<usize>,
    dummy: ::std::marker::PhantomData<G>
}

impl<G> GateBuilder<G> where G: ?Sized
{
    fn new(param_idxs: &[usize]) -> Self
    {
        GateBuilder
        {
            param_idxs: param_idxs.to_owned(),
            dummy: ::std::marker::PhantomData
        }
    }

    fn build(&self, params: &[f64]) -> G
    where G: gates::Gate
    {
        let n = self.param_idxs.len();
        match n
        {
            0 => G::new_parametrized(&[]),
            1 => G::new_parametrized(&[params[self.param_idxs[0]]]),
            2 => G::new_parametrized(&[params[self.param_idxs[0]], params[self.param_idxs[1]]]),
            _ => {
                let mut v = vec![];
                for &idx in self.param_idxs.iter()
                {
                    v.push(params[idx]);
                }
                G::new_parametrized(&v)
            }
        }
    }
}

/// Operation in a Custom gate.
enum CustomDefOp
{
    /// Unary gate operating on a single qubit.
    Unary(Box<GateBuilder<gates::UnaryGate>>, usize),
    /// Binary gate operating on two qubits.
    Binary(GateBuilder<gates::BinaryGate>, usize, usize),
    /// Gate operating on multiple qubits.
    Nary(GateBuilder<gates::NaryGate>, Vec<usize>)
}

impl CustomDef
{
    /// Create a new custom gate.
    ///
    /// Initialize a new custom gate definition with name `name` for operating
    /// on `nr_bits` qubits at a time. The gate types making up the operation
    /// should be added using the `add_unary_gate()`, `add_binary_gate()` or
    /// `add_n_ary_gate()` functions.
    pub fn new(name: &str, nr_bits: usize) -> Self
    {
        CustomDef
        {
            name: name.to_owned(),
            nr_bits: nr_bits,
            nr_params: 0,
            ops: vec![]
        }
    }

    /// Add a unary gate.
    ///
    /// Append an unparametrized unary gate `gate` operating on qubit `bit` to
    // this gate definition.
    pub fn add_unary_gate<G: 'static>(&mut self, bit: usize)
    where G: gates::UnaryGate
    {
        assert!(bit < self.nr_bits, "Invalid bit index {} for {}-bit gate", bit, self.nr_bits);
        self.ops.push(CustomDefOp::Unary(Box::new(GateBuilder::new(&[])), bit));
    }

    /// Add a binary gate.
    ///
    /// Append an  unparamaterixed binary gate `gate` operating on qubits
    /// `bit0` and `bit1` to this gate.
    pub fn add_binary_gate<G: 'static>(&mut self, gate: G, bit0: usize, bit1: usize)
    where G: gates::BinaryGate
    {
        assert!(bit0 < self.nr_bits, "Invalid bit index {} for {}-bit gate", bit0, self.nr_bits);
        assert!(bit1 < self.nr_bits, "Invalid bit index {} for {}-bit gate", bit1, self.nr_bits);
        self.ops.push(CustomDefOp::Binary(GateBuilder::new(&[]), bit0, bit1));
    }

    /// Add a multi-bit gate.
    ///
    /// Append an unparametrized `n`-ary gate `gate`, operating on the `n`
    /// qubits in `bits`, to this gate.
    pub fn add_n_ary_gate<G: 'static>(&mut self, gate: G, bits: &[usize])
    where G: gates::NaryGate
    {
        for &bit in bits.iter()
        {
            assert!(bit < self.nr_bits, "Invalid bit index {} for {}-bit gate", bit, self.nr_bits);
        }
        self.ops.push(CustomDefOp::Nary(GateBuilder::new(&[]), bits.to_owned()));
    }

    pub fn build_gate(&self, params: &[f64]) -> gates::Custom
    {
        let gate = gates::Custom::new(&self.name, self.nr_bits);
        for op in self.ops.iter()
        {
            match *op
            {
                CustomDefOp::Unary(ref builder, bit) => {
                    let g = builder.build(&[]);
                    gate.add_unary_gate(g, bit);
                }
                _ => { }
            }
        }
        gate
    }
}

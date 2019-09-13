/// Trait for a quantum state
///
/// Trait QuState describes the interface for changing, measuring and monitoring
/// the quantum state of the simulated computer.
pub trait QuState
{
    /// Apply a n-ary quantum gate `gate` on the qubits from `bits` in this state.
    fn apply_gate<G>(&mut self, gate: &G, bits: &[usize]) -> crate::error::Result<()>
    where G: crate::gates::Gate + ?Sized;
    /// Apply a unary gate to all qubits
    ///
    /// Apply the single-bit gate `gate` to all qubits in the quantum state.
    fn apply_unary_gate_all<G>(&mut self, gate: &G) -> crate::error::Result<()>
    where G: crate::gates::Gate + ?Sized;
    /// Apply a conditional n-ary quantum gate `gate`, controlled by classical
    /// bit `control`, on the qubits from `bits` in this state.
    fn apply_conditional_gate<G>(&mut self, control: &[bool], gate: &G,
        bits: &[usize]) -> crate::error::Result<()>
    where G: crate::gates::Gate + ?Sized;

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `qbit` in the state. Measurement is done
    /// in the `z`-basis. The random number generator `rng` is used for sampling.
    /// The result is returned as an array containing the measurement result for
    /// each run.
    fn measure<R: rand::Rng>(&mut self, qbit: usize, rng: &mut R)
        -> crate::error::Result<ndarray::Array1<u64>>;
    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `qbit` in the state to classical bit
    /// `cbit` in `res`, which should be an array of sufficient length to store
    /// results for the total number of runs in the state. Measurement is done
    /// in the `z`-basis. The random number generator `rng` is used for sampling.
    fn measure_into<R: rand::Rng>(&mut self, qbit: usize, cbit: usize,
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>;
    /// Measure all qubits
    ///
    /// Measure all qubits in this state, and return the results. The random
    /// number generator `rng` is used for sampling.
    fn measure_all<R: rand::Rng>(&mut self, rng: &mut R)
        -> crate::error::Result<ndarray::Array1<u64>>;
    /// Measure all qubits
    ///
    /// Measure all qubits in this state, and store the results in `res`,
    /// which should be of sufficient length to hold results for the number of
    /// runs in this state.  The first qubit measured is stored at the bit
    /// position indicated by the first element of `cbits`, and so on. The random
    /// number generator `rng` is used for sampling.
    fn measure_all_into<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>;

    /// Measure a qubit.
    ///
    /// Perform a measurement on qubit `qbit` in the state to bit `cbit` in res,
    /// without affecting the quantum state, i.e. the wave function is not
    /// collapsed. The output array `res` should be of sufficient length to store
    /// results for the total number of runs in the state. Measurement is done
    /// in the `z`-basis. The random number generator `rng` is used for sampling.
    /// NOTE: this is not a physical process, and impossible to reproduce on
    /// a real quantum computer.
    fn peek_into<R: rand::Rng>(&self, qbit: usize, cbit: usize,
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>;
    /// Measure all qubits
    ///
    /// Measure all qubits in this state without affecting the quantum state,
    /// i.e. without collapsing the wave function. The measurement results
    /// are stored in `res`, which must be of sufficient length to hold results
    /// for the total number of runs in the state. The first qubit measured
    /// is stored at the bit position indicated by the first element of `cbits`,
    /// and so on. The random number generator `rng` is used for sampling.
    /// NOTE: this is not a physical process, and impossible to reproduce on
    /// a real quantum computer.
    fn peek_all_into<R: rand::Rng>(&mut self, cbits: &[usize],
        res: &mut ndarray::Array1<u64>, rng: &mut R) -> crate::error::Result<()>;

    /// Reset a qubit
    ///
    /// Reset the qubit with index `bit` to zero. This is done by measuring the
    /// bit, and rotating it back to zero if the result is 1. The random
    /// number generator `rng` is used for sampling in the measurement.
    fn reset<R: rand::Rng>(&mut self, bit: usize, rng: &mut R)
        -> crate::error::Result<()>;
    /// Reset all qubits
    ///
    /// Reset all qubits in this experiment, returning the state to |00...0⟩
    /// for all runs.
    fn reset_all(&mut self);
}

/// Collect which states to apply conditional gate to into ranges
///
/// A quantum calculation is represented by a set of quantum states, where each
/// separate state is combined with the number of times it has been attained
/// in the calculation. Given such a quantum state, and `N` control bits
/// (where `N` is the total number of runs in the calculation), this function
/// collects triplets (`start`, `length`, `apply`), where a conditional
/// gate should be applied to `length` copies of the quantum state at index
/// `start` if `apply` is `true`, and not copied unchanged if it is `false`.
pub fn collect_conditional_ranges(counts: &[usize], control: &[bool])
    -> Vec<(usize, usize, bool)>
{
    let mut ranges = vec![];
    let mut off = 0;
    for (icol, &count) in counts.iter().enumerate()
    {
        let mut begin = off;
        let mut prev = control[off];
        for ibit in off+1..off+count
        {
            if control[ibit] != prev
            {
                ranges.push((icol, ibit-begin, prev));
                begin = ibit;
                prev = !prev;
            }
        }
        if begin < off+count
        {
            ranges.push((icol, off+count-begin, prev));
        }

        off += count;
    }

    ranges
}

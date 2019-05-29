# [0.4.0] 2019-05-29

Much more efficient simulation of stabilizer circuits. Stabilizer circuits consist
of Clifford gates only: gates that can be described entirely in terms of `H`, `S`
and `CX` gates. These circuits can be represented more efficiently, scaling
quadratically in the number of qubits instead of exponentially as when using the
standard coefficient vector representation (see [ArXiv:9807006](https://arxiv.org/abs/quant-ph/9807006)
for details). The gates included with `q1tsim` that can be simulated in this
representation are `H`, `X`, `Y`, `Z`, `S` and `S`<sup>`†`</sup>,
`V` and `V`<sup>`†`</sup>, `CX`, `CY`, `CZ`, and `Swap`,
as well as `Composite` gates containing only these gates.  With the new backend,
simulation of stabilizer circuits with thousands of qubits is possible.

By default, `q1tsim` will now use the stabilizer representation if possible. If
you wish to use the coefficient vector representation instead when simulating a
stabilizer circuit, use the `Circuit::execute_with()` method.

# [0.3.0] 2019-04-09

- Update code to use Rust 2018 edition.
- Most methods of `Circuit` now return a `Result` with a possible error code.
Checking the result of every call while building a circuit may be tedious. In that
case, use the new `circuit!` macro to build the circuit, which stops at the first
error it encounters, e.g.
```
let circuit = circuit!(2, 2, {
    h(0);
    x(1);
    measure_all(&[0, 1]);
}).expect("Failed to build circuit");
```
- The number of runs to average over is no longer specified in `Circuit::new()`,
but passed in later when calling `Circuit::execute()`.
- A new `Circuit::reexecute()` method was added, to execute an existing circuit
again without resetting the state, i.e. starting in the end state of the previous
execution.
- The `Gate` trait no longer requires you to implement the `cost()` function.
Its application is limited, and sometimes hard to estimate. The function now
returns infinity by default.
- New `peek()` and `peek_all()` methods on `Circuit` that allow you to perform
a measurement without collapsing the quantum state. This is of course not a
physical process, and cannot be reproduced on a quantum computer, but can be
useful for simualting intermediate measurements while still evolving the
quantum state.
- The random number generator used for sampling in measurements is no longer
stored in `QuState`, but instead passed in. Combined with the new `execute_with_rng()`
and reexecute_with_rng()` methods on `Circuit`, this allows users to use their
own RNG for sampling. This can be useful, e.g., to create reproducible simulations.
- The internal representation of QuState was simplified to use a single matrix
of quantum states instead of a loose collection of individual states. This allows
us to perform operations on all states at once. Calculations with heavy branching
(due to measurements), should run much faster now.
- Several updates and fixes in export to various formats.
- New `s()` and `sdg` methods on `Circuit` for quickly adding `S` or `S`<sup>`†`</sup>
phase gates to a circuit.
- Added doubly controlled rotation gates `CCRX`, `CCRY`, and `CCRZ`.


# [0.2.0] 2019-02-27

- Much better test coverage
- Export of circuits to `LaTeX` code, using the `qcircuit` package.
- Arithmetic expressions, and the constant `pi` can now be used in describing
gate parameters in `Composite::from_string()`.
- Measurements off all bits can now be performed in a different basis than just
the Pauli `Z` basis.
- Several small performance improvements
- Fixed wrong bit order of the control register in classically controlled
operations. The bit order is now consistent with the histograms from
measurements (i.e. first control bit is now least significant).
- Several bug fixes in the export to OpenQasm en c-Qasm.
- Fixed parameters being passed in the wrong order by `Circuit::u3()`.
- Removed `Gate::expanded_matrix()`. It was used nowhere, and it can easily be
reimplemented more generally by operating with the gate on a matrix. Use the
identity matrix to obtain the effect of the original function.

# [0.1.0] 2019-01-30

Initial version

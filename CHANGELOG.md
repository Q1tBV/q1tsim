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

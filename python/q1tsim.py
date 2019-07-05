import cffi
import json
import q1tsimffi

class RefParam(object):
    """Structure for reference parameters

    A RefPrama struct is used for passing parameters to a gate by reference,
    i.e. one can change these parameters between executions of the same circuit.
    To use them, simply construct a RefParam with an initial value

    theta = RefParam(3.14)

    and pass it to a gate expecting a float parameter in the same way a direct
    value would be given

    circuit.cx(theta)

    To change the parameter, use the assign() method. The next execution of the
    circuit will use the updated value, without constructing a new circuit.
    """
    def __init__(self, value):
        """Create a new reference parameter with initial value value."""
        self.__ptr = cffi.FFI().new('double *', value)

    def __float__(self):
        """Convert to float, i.e. return the value of this parameter."""
        return self.__ptr[0]

    def assign(self, value):
        """Asign a new value value to this parameter."""
        self.__ptr[0] = value

    def pointer(self):
        """Return the pointer to the parameter value."""
        return self.__ptr

class Circuit(object):
    """A quantum circuit

    Struct Circuit represents a quantum circuit, holding a quantum state and the
    operations to be performed on it.
    """

    def __init__(self, nr_qbits, nr_cbits=0):
        """Create a new circuit.

        Create a new (empty) quantum circuit, with nr_qbits quantum bits and
        nr_cbits classical bits.
        """
        self.__sim = q1tsimffi.q1tsim()
        self.__ptr = self.__sim.circuit_new(nr_qbits, nr_cbits)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__sim.circuit_free(self.__ptr)
        self.__ptr = None

    def nr_qbits(self):
        """The number of quantum bits in this circuit"""
        return self.__sim.circuit_nr_qbits(self.__ptr)

    def nr_cbits(self):
        """The number of classical bits in this circuit"""
        return self.__sim.circuit_nr_cbits(self.__ptr)

    def add_gate(self, name, qbits, params=None):
        """Add a gate.

        Append a n-ary quantum gate gate, operating on the n qubits in bits, to
        this circuit.
        """
        cname = bytes(name, 'utf-8')
        if params is None:
            res = q1tsimffi.unpack_result(
                self.__sim.circuit_add_gate(self.__ptr, cname, qbits, len(qbits), cffi.FFI.NULL, 0)
            )
        else:
            cparams = q1tsimffi.make_parameters(params)
            res = q1tsimffi.unpack_result(
                self.__sim.circuit_add_gate(self.__ptr, cname, qbits, len(qbits), cparams, len(params))
            )
        return res

    def add_conditional_gate(self, control, target, name, qbits, params=None):
        """Add a conditional gate.

        Append a n-ary gate gate, that will operate on the n qubits in
        bits to this circuit. The gate will only be applied only when the
        classical bits with indices from control form the target word target.
        The bit at the position of the first index in control is interpreted
        as the least significant bit to check.
        """
        cname = bytes(name, 'utf-8')
        if params is None:
            res = q1tsimffi.unpack_result(
                self.__sim.circuit_add_conditional_gate(self.__ptr,
                    control, len(control), target, cname, qbits, len(qbits),
                    cffi.FFI.NULL, 0
                )
            )
        else:
            res = q1tsimffi.unpack_result(
                self.__sim.circuit_add_conditional_gate(self.__ptr,
                    control, len(control), target, cname, qbits, len(qbits),
                    params, len(params)
                )
            )
        return res

    def ch(self, control, target):
        """Add a controlled Hadamard gate.

        Add a controlled Hadamard gate, controlled by qubit control, and
        operating on qubit target, to this circuit.
        """
        return self.add_gate('CH', [control, target])

    def crx(self, theta, control, target):
        """Add a conditional RX gate.

        Add a conditional RX(θ) gate, controlled by qubit control, and operating
        on qubit target, to this circuit.
        """
        return self.add_gate('CRX', [control, target], [theta])

    def cry(self, theta, control, target):
        """Add a conditional RY gate.

        Add a conditional RY(θ) gate, controlled by qubit control, and operating
        on qubit target, to this circuit.
        """
        return self.add_gate('CRY', [control, target], [theta])

    def crz(self, lmb, control, target):
        """Add a conditional RZ gate.

        Add a conditional RZ(λ) gate, controlled by qubit control, and operating
        on qubit target, to this circuit.
        """
        return self.add_gate('CRZ', [control, target], [theta])

    def cx(self, control, target):
        """Add a CX gate.

        Add a controlled X gate (controlled NOT), controlled by qubit control,
        and operating on qubit target, to this circuit.
        """
        return self.add_gate('CX', [control, target])

    def cy(self, control, target):
        """Add a CY gate.

        Add a controlled Y gate, controlled by qubit control, and
        operating on qubit target, to this circuit.
        """
        return self.add_gate('CY', [control, target])

    def cz(self, control, target):
        """Add a CZ gate.

        Add a controlled Z gate, controlled by qubit control, and
        operating on qubit target, to this circuit.
        """
        return self.add_gate('CZ', [control, target])

    def h(self, qbit):
        """Add a Hadamard gate.

        Add a Hadamard gate operating on qubit `qbit`, to this circuit.
        """
        return self.add_gate('H', [qbit])

    def i(self, qbit):
        """Add an identity gate.

        Add an identity gate operating on qubit `qbit`, to this circuit. Since
        this gate does nothing, you might want to consider if you really need
        it.
        """
        return self.add_gate('I', [qbit])

    def rx(self, theta, qbit):
        """Add a RX gate.

        Add an RX(θ) gate operating on qubit `bit`, to this circuit.
        """
        return self.add_gate('RX', [qbit], [theta])

    def ry(self, theta, qbit):
        """Add a RY gate.

        Add an RY(θ) gate operating on qubit `bit`, to this circuit.
        """
        return self.add_gate('RY', [qbit], [theta])

    def rz(self, lmb, qbit):
        """Add a RZ gate.

        Add an RZ(λ) gate operating on qubit `bit`, to this circuit.
        """
        return self.add_gate('RZ', [qbit], [lmb])

    def s(self, qbit):
        """Add a phase gate

        Add an S phase gate (rotation of π/2 around the Z axis) operating on
        qubit bit, to this circuit.
        """
        return self.add_gate('S', [qbit])

    def sdg(self, qbit):
        """Add an inverse phase gate

        Add an S† gate, the inverse of the S gate, operating on qubit bit,
        to this circuit.
        """
        return self.add_gate('Sdg', [qbit])

    def swap(self, qbit0, qbit1):
        """Add a swap gate.

        Add a swap gate, swapping qubits qbit0 and qbit1.
        """
        return self.add_gate('Swap', [qbit0, qbit1])

    def t(self, qbit):
        """Add a T gate

        Add an T phase gate (rotation of π/4 around the Z axis) operating on
        qubit bit, to this circuit.
        """
        return self.add_gate('T', [qbit])

    def tdg(self, qbit):
        """Add an inverse T gate

        Add an T† gate, the inverse of the T gate, operating on qubit bit,
        to this circuit.
        """
        return self.add_gate('Tdg', [qbit])

    def u1(self, lmb, qbit):
        """Add a U1 gate.

        Add a U1(λ) gate operating on qubit qbit, to this circuit. This gate is,
        up to a global phase, equivalent to the RZ gate.
        """
        return self.add_gate('U1', [qbit], [lmb])

    def u2(self, phi, lmb, qbit):
        """Add a U2 gate.

        Add a U2(ϕ, λ) gate operating on qubit qbit, to this circuit.
        """
        return self.add_gate('U2', [qbit], [phi, lmb])

    def u3(self, theta, phi, lmb, qbit):
        """Add a U3 gate.

        Add a U3(θ, ϕ, λ) gate operating on qubit qbit, to this circuit.
        """
        return self.add_gate('U3', [qbit], [theta, phi, lmb])

    def v(self, qbit):
        """Add a V gate.

        Add a V gate (square root of NOT) operating on qubit qbit, to this circuit.
        """
        return self.add_gate('V', [qbit])

    def vdg(self, qbit):
        """Add an inverse V gate.

        Add an V† gate, the inverse of the V gate, operating on qubit bit,
        to this circuit.
        """
        return self.add_gate('Vdg', [qbit])

    def x(self, qbit):
        """Add a Pauli X gate.

        Add a Pauli X gate (NOT gate) operating on qubit qbit, to this circuit.
        """
        return self.add_gate('X', [qbit])

    def y(self, qbit):
        """Add a Pauli Y gate.

        Add a Pauli Y gate operating on qubit qbit, to this circuit.
        """
        return self.add_gate('Y', [qbit])

    def z(self, qbit):
        """Add a Pauli Z gate.

        Add a Pauli Z gate operating on qubit qbit, to this circuit.
        """
        return self.add_gate('Z', [qbit])

    def __peek_measure_basis(self, qbit, cbit, basis, collapse):
        cbasis = bytes(basis, 'utf-8')
        ccollapse = 1 if collapse else 0
        res = q1tsimffi.unpack_result(
            self.__sim.circuit_measure(self.__ptr, qbit, cbit, cbasis, ccollapse)
        )
        return res

    def __peek_measure_all_basis(self, cbits, basis, collapse):
        cbasis = bytes(basis, 'utf-8')
        ccollapse = 1 if collapse else 0
        res = q1tsimffi.unpack_result(
            self.__sim.circuit_measure_all(self.__ptr, cbits, len(cbits), cbasis, ccollapse)
        )
        return res

    def measure_basis(self, qbit, cbit, basis):
        """Add a measurement

        Add measurement of qubit qbit in basis basis, into classical bit
        cbit, to this circuit.
        """
        return self.__peek_measure_basis(qbit, cbit, basis, True)

    def measure_x(self, qbit, cbit):
        """Add a measurement.

        Add measurement of qubit qbit in the Pauli X basis, into classical
        bit cbit to this circuit.
        """
        return self.measure_basis(qbit, cbit, 'X')

    def measure_y(self, qbit, cbit):
        """Add a measurement.

        Add measurement of qubit qbit in the Pauli Y basis, into classical
        bit cbit to this circuit.
        """
        return self.measure_basis(qbit, cbit, 'Y')

    def measure_z(self, qbit, cbit):
        """Add a measurement.

        Add measurement of qubit qbit in the Pauli Z basis, into classical
        bit cbit to this circuit.
        """
        return self.measure_basis(qbit, cbit, 'Z')

    def measure(self, qbit, cbit):
        """Add a measurement.

        Add measurement of qubit qbit into classical bit cbit to this circuit.
        This is an alias for measure_z().
        """
        return self.measure_z(qbit, cbit)

    def measure_all_basis(self, cbits, basis):
        """Add a measurement.

        Add the measurement of all qubits in the quantum state into the classical
        bits cbits. Measurement is done in basis basis.
        """
        return self.__peek_measure_all_basis(cbits, basis, True)

    def measure_all(self, cbits):
        """Add a measurement.

        Add the measurement of all qubits in the quantum state into the classical
        bits cbits. Measurement is done in the Pauli Z basis.
        """
        return self.measure_all_basis(cbits, 'Z')

    def peek_basis(self, qbit, cbit, basis):
        """Add a measurement.

        Add the measurement of qubit qbit in the quantum state into the
        classical bit cbit. Measurement is done in basis basis, without
        collapsing the quantum state.
        NOTE: this is not a physical process, and cannot be reproduced on a real
        quantum computer.
        """
        return self.__peek_measure_basis(qbit, cbit, basis, False)

    def peek_x(self, qbit, cbit):
        """Add a measurement.

        Add the measurement of qubit qbit in the quantum state into the
        classical bit cbit. Measurement is done in the Pauli X basis, without
        collapsing the quantum state.
        NOTE: this is not a physical process, and cannot be reproduced on a real
        quantum computer.
        """
        return self.peek_basis(qbit, cbit, 'X')

    def peek_y(self, qbit, cbit):
        """Add a measurement.

        Add the measurement of qubit qbit in the quantum state into the
        classical bit cbit. Measurement is done in the Pauli Y basis, without
        collapsing the quantum state.
        NOTE: this is not a physical process, and cannot be reproduced on a real
        quantum computer.
        """
        return self.peek_basis(qbit, cbit, 'Y')

    def peek_z(self, qbit, cbit):
        """Add a measurement.

        Add the measurement of qubit qbit in the quantum state into the
        classical bit cbit. Measurement is done in the Pauli Z basis, without
        collapsing the quantum state.
        NOTE: this is not a physical process, and cannot be reproduced on a real
        quantum computer.
        """
        return self.peek_basis(qbit, cbit, 'Z')

    def peek(self, qbit, cbit):
        """Add a measurement.

        Add the measurement of qubit qbit in the quantum state into the
        classical bit cbit. This is an alias foor peek_z().
        """
        return self.peek_z(qbit, cbit)

    def peek_all_basis(self, cbits, basis):
        """Add a measurement.

        Add the measurement of all qubits in the quantum state into the classical
        bits cbits. Measurement is done in basis basis, without
        collapsing the quantum state.
        NOTE: this is not a physical process, and cannot be reproduced on a real
        quantum computer.
        """
        return self.__peek_measure_all_basis(cbits, basis, False)

    def peek_all(self, cbits):
        """Add a measurement.

        Add the measurement of all qubits in the quantum state into the classical
        bits cbits. Measurement is done in the Pauli `Z` basis, without
        collapsing the quantum state.
        NOTE: this is not a physical process, and cannot be reproduced on a real
        quantum computer.
        """
        return self.peek_all_basis(cbits, 'Z')

    def reset(self, qbit):
        """Reset a qubit

        Reset the qubit qbit to |0⟩. This is done by measuring the bit, and
        flipping it if the result is 1, so this is potentially an expensive
        operation.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_reset(self.__ptr, qbit))
        return res

    def reset_all(self):
        """Reset all qubits.

        Reset the entire quantum state of the circuit to |00...0⟩. The classical
        register is not affected.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_reset_all(self.__ptr))
        return res

    def execute(self, nr_shots):
        """Execute this circuit

        Execute this circuit, performing its operations and measurements.
        Measurements are made over nr_shots executions of the circuit. This
        function clears any previous states of the system (quantum or classical).
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_execute(self.__ptr, nr_shots))
        return res

    def reexecute(self):
        """Execute a circuit again.

        Run this circuit again, starting with the state from the previous
        execution.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_reexecute(self.__ptr))
        return res

    def histogram(self):
        """Create a histogram of measurements.

        Create a histogram of the measured classical bits and return it as a
        dictionary mapping measurement result to the number of times the result
        was measured. The n bits in the classical register are collected in a
        string key, with the last character in the key corresponding to the
        first bit (at index 0) in the classical register and vice versa.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_histogram(self.__ptr))
        return res

    def latex(self):
        """Export to OpenQasm

        Export this circuit to LaTeX using the qcircuit pacakge. On a successful
        conversion, returns a string with the LaTeX code.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_latex(self.__ptr))
        return res

    def open_qasm(self):
        """Export to OpenQasm

        Export this circuit to a program in OpenQasm format. On a successful
        conversion, returns a string with the program text.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_open_qasm(self.__ptr))
        return res

    def c_qasm(self):
        """Export to c-Qasm.

        Export this circuit to a program in c-Qasm format. On a successful
        conversion, returns a string with the program text.
        """
        res = q1tsimffi.unpack_result(self.__sim.circuit_c_qasm(self.__ptr))
        return res

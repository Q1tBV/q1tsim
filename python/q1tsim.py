import cffi
import json
import q1tsimffi

class Circuit(object):
    def __init__(self, nr_qbits, nr_cbits=0):
        self.__sim = q1tsimffi.q1tsim()
        self.__ptr = self.__sim.circuit_new(nr_qbits, nr_cbits)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__sim.circuit_free(self.__ptr)
        self.__ptr = None

    def nr_qbits(self):
        return self.__sim.circuit_nr_qbits(self.__ptr)

    def nr_cbits(self):
        return self.__sim.circuit_nr_cbits(self.__ptr)

    def add_gate(self, name, qbits, params=None):
        cname = bytes(name, 'utf-8')
        if params is None:
            res = q1tsimffi.unpack_result(
                self.__sim.circuit_add_gate(self.__ptr, cname, qbits, len(qbits), cffi.FFI.NULL, 0)
            )
        else:
            res = q1tsimffi.unpack_result(
                self.__sim.circuit_add_gate(self.__ptr, cname, qbits, len(qbits), params, len(params))
            )
        return res

    def ch(self, control, target):
        return self.add_gate('CH', [control, target])

    def cx(self, control, target):
        return self.add_gate('CX', [control, target])

    def cy(self, control, target):
        return self.add_gate('CY', [control, target])

    def cz(self, control, target):
        return self.add_gate('CZ', [control, target])

    def h(self, qbit):
        return self.add_gate('H', [qbit])

    def i(self, qbit):
        return self.add_gate('I', [qbit])

    def rx(self, theta, qbit):
        return self.add_gate('RX', [qbit], [theta])

    def ry(self, theta, qbit):
        return self.add_gate('RY', [qbit], [theta])

    def rz(self, lmb, qbit):
        return self.add_gate('RZ', [qbit], [theta])

    def s(self, qbit):
        return self.add_gate('S', [qbit])

    def sdg(self, qbit):
        return self.add_gate('Sdg', [qbit])

    def swap(self, qbit0, qbit1):
        return self.add_gate('Swap', [qbit0, qbit1])

    def t(self, qbit):
        return self.add_gate('T', [qbit])

    def tdg(self, qbit):
        return self.add_gate('Tdg', [qbit])

    def u1(self, lmb, qbit):
        return self.add_gate('U1', [qbit], [lmb])

    def u2(self, phi, lmb, qbit):
        return self.add_gate('U2', [qbit], [phi, lmb])

    def u3(self, theta, phi, lmb, qbit):
        return self.add_gate('U3', [qbit], [theta, phi, lmb])

    def v(self, qbit):
        return self.add_gate('V', [qbit])

    def vdg(self, qbit):
        return self.add_gate('Vdg', [qbit])

    def x(self, qbit):
        return self.add_gate('X', [qbit])

    def y(self, qbit):
        return self.add_gate('Y', [qbit])

    def z(self, qbit):
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
        self.__peek_measure_basis(qbit, cbit, basis, True)

    def measure_x(self, qbit, cbit):
        return self.measure_basis(qbit, cbit, 'X')

    def measure_y(self, qbit, cbit):
        return self.measure_basis(qbit, cbit, 'Y')

    def measure_z(self, qbit, cbit):
        return self.measure_basis(qbit, cbit, 'Z')

    def measure(self, qbit, cbit):
        return self.measure_z(qbit, cbit)

    def measure_all_basis(self, cbits, basis):
        self.__peek_measure_all_basis(cbits, basis, True)

    def measure_all(self, cbits):
        return self.measure_all_basis(cbits, 'Z')

    def peek_basis(self, qbit, cbit, basis):
        self.__peek_measure_basis(qbit, cbit, basis, False)

    def peek_x(self, qbit, cbit):
        return self.peek_basis(qbit, cbit, 'X')

    def peek_y(self, qbit, cbit):
        return self.peek_basis(qbit, cbit, 'Y')

    def peek_z(self, qbit, cbit):
        return self.peek_basis(qbit, cbit, 'Z')

    def peek(self, qbit, cbit):
        return self.peek_z(qbit, cbit)

    def peek_all_basis(self, cbits, basis):
        self.__peek_measure_all_basis(cbits, basis, False)

    def peek_all(self, cbits):
        return self.peek_all_basis(cbits, 'Z')

    def execute(self, nr_shots):
        res = q1tsimffi.unpack_result(self.__sim.circuit_execute(self.__ptr, nr_shots))
        return res

    def histogram(self):
        res = q1tsimffi.unpack_result(self.__sim.circuit_histogram(self.__ptr))
        return res

    def latex(self):
        res = q1tsimffi.unpack_result(self.__sim.circuit_latex(self.__ptr))
        return res

    def open_qasm(self):
        res = q1tsimffi.unpack_result(self.__sim.circuit_open_qasm(self.__ptr))
        return res

    def c_qasm(self):
        res = q1tsimffi.unpack_result(self.__sim.circuit_c_qasm(self.__ptr))
        return res

import q1tsim

nr_qbits = 2
nr_cbits = 2
nr_shots = 1000

with q1tsim.Circuit(nr_qbits, nr_cbits) as circuit:
    circuit.h(0)
    circuit.h(1)
    circuit.measure(0, 0)
    circuit.measure(1, 1)
    circuit.execute(nr_shots)
    print(circuit.cstate())

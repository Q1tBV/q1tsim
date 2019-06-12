import q1tsim

nr_qbits = 1
nr_cbits = 1
nr_shots = 10**3

ang = 1.23
with q1tsim.Circuit(nr_qbits, nr_cbits) as circuit:
    circuit.rx(ang, 0)
    circuit.peek(0, 0)
    circuit.execute(nr_shots)
    print(circuit.histogram())

    ang = -1.23
    circuit.reexecute()
    print(circuit.histogram())

print()

ang = q1tsim.RefParam(1.23)
with q1tsim.Circuit(nr_qbits, nr_cbits) as circuit:
    circuit.rx(ang, 0)
    circuit.peek(0, 0)
    circuit.execute(nr_shots)
    print(circuit.histogram())

    ang.assign(-1.23)
    circuit.reexecute()
    print(circuit.histogram())

import os
import sys
import numpy as np

for line in sys.stdin.read().splitlines():
    line = line.strip()
    if line == "" or line.startswith(("//", "OPENQASM", "include")):
        continue
    elif line.startswith("qreg"):
        num_qubits = int(line.split("[")[1].split("]")[0])
        parity = np.identity(num_qubits, dtype=int)
    elif line.startswith("cx"):
        parts = line.split("[")
        q1 = int(parts[1].split("]")[0])
        q2 = int(parts[2].split("]")[0])
        parity[q2] ^= parity[q1]
    else:
        print("QASM line not recognized:", line)
        exit(-1)

with open("parity_matrix.txt", "w") as pmfile:
    for x in parity:
        for y in x:
            print(y, end=" ", file=pmfile)
        print(file=pmfile)

os.system(f"sh matrix_cnot.sh -Q{num_qubits} parity_matrix.txt")

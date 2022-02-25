import numpy as np
import math
from qiskit import QuantumCircuit

def qft_dagger(qc, n):
    """
    extracted from opensource qiskit textbook
    https://qiskit.org/textbook/ch-algorithms/quantum-phase-estimation.html

    n-qubit QFTdagger the first n qubits in circ"""
    # Don't forget the Swaps!
    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        for m in range(j):
            qc.cp(-math.pi/float(2**(j-m)), m, j)
        qc.h(j)
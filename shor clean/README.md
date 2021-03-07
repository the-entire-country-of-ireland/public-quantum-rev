This is a unobfuscated implementation of Shor's algorithm for factoring integers in qiskit.
Uses `2*n + n_e + 3` qubits.

Intended for educational purposes; mainly to visualize the output distribution. Try running the provided script, varying `ne` from 6 to 10, and watch how the peaks get sharper and sharper as you increase the precision of the quantum phase estimation algorithm.

Uses a QFT-based adder, which substantially reduces the required number of ancilla qubits.
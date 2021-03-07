import numpy as np
from fractions import Fraction
from matplotlib import pyplot as plt

from math import gcd, log2, pi
import random

# Import Qiskit
from qiskit import QuantumCircuit
from qiskit import Aer, execute


########################################

# circuit parameters

"""
to match Figure 5.1 in original Shor paper, use:
N = 33
a = 5
n = 6
ne = 8

another good example is:
N = 55
a = 3
n = 7
ne = 8
"""

# number to factor
N = 33

# base of exponent
a = 5

print("a = {}\nN = {}\n".format(a,N))
assert gcd(a,N) == 1
assert 2 <= a < N


########################################

# number of bits + padding to store N
# you need an extra bit or 2 of padding to prevent overflow.
# 2 bits of padding always works, 1 sometimes works
n = 6

# number of bits in exponent
# = number of classical bits as output.
# the more exponent bits, the sharper the peaks in the output
# --> the more precise the estimate of the order.
ne = 6


assert n - 0.5 > log2(N), "increase the number of n bits"


# total number of qubits used
m = 2*n + ne + 1

# stores result
x_qubits = list(range(ne, n+ne))
# auxillary bits used for controlled multiplication
y_qubits = list(range(n+ne, 2*n + ne))
# exponent bits
z_qubits = list(range(ne))
# auxillary bit used in modular addition to copy the carry/overflow bit
ancillary_qubit = m - 1


########################################


# we represent negative numbers via the ones-complement

# little-endian encoding (LSB first)
# ones complement
def int_to_bin_LSB(x,n):
    sign = (x < 0)
    x = abs(x)
    x = bin(x)[2:]
    x = x.zfill(n)[:n]

    x = x[::-1]
    x = [int(i) for i in x]
    return x


# LSB first
def bin_to_int_LSB(L):
    x = ''.join([str(i) for i in L])
    x = eval('0b' + x[::-1])
    return x


def ones_complement(x,n):
    return (1<<n) + ~x


########################################



# classical constant factors, in little endian (LSB first)
def get_const_factors(a,N,n):
    const_factors = [(a<<i)%N for i in range(n)]
    const_factors = [int_to_bin_LSB(t,n) for t in const_factors]
    const_factors = [[int(i) for i in t] for t in const_factors]

    return const_factors


def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def mod_inverse(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

########################################

def cphi_gate(circ, k, control, target):
    theta = 2 * pi / (1<<k)
    circ.cp(theta, control, target)

def phi_gate(circ, k, target):
    theta = 2 * pi / (1<<k)
    circ.p(theta, target)


# negation circuit
def ones_complement_circuit(a_qubits, m):
    circ = QuantumCircuit(m,ne)
    
    for i in a_qubits:
        circ.x(i)
        
    return circ

    
########################################
    
# initialization circuit

def make_initialization_circuit(a_qubits, a_bits, m):
    initialize_circ = QuantumCircuit(m,ne)
    for i,v in enumerate(a_bits):
        if v:
            initialize_circ.x(a_qubits[i])
            
    return initialize_circ


########################################
    
# quantum fourier transform
        
QFT_cutoff = 5 + int(log2(n))

def make_QFT_circ(a_qubits,n,m):
    
    QFT_circ = QuantumCircuit(m,ne)
    
    QFT_cutoff = 5 + int(log2(n))
    
    for i in range(n)[::-1]:
        QFT_circ.h(a_qubits[i])
        
        for j in range(i)[::-1]:
            k = i - j + 1
            
            if k > QFT_cutoff:
                break
            
            cphi_gate(QFT_circ, k,a_qubits[j],a_qubits[i])
            
    return QFT_circ



########################################

# quantum adder
# in-place a+b, where a is quantum and b is classical

def make_adder_circ(b_bits, a_qubits, m):
    adder_circ = QuantumCircuit(m,ne)
    
    for i in range(n)[::-1]:
        
        theta = 0
        
        for j in range(i + 1)[::-1]:
            k = i - j + 1
            
            if k > QFT_cutoff:
                break
            
            if b_bits[j]:
                theta += 2 * pi / (1<<k)
        
        if theta != 0:
            adder_circ.p(theta, a_qubits[i])
            
    return adder_circ

def make_controlled_adder_circ(control_qubit, b_bits, a_qubits, m):
    adder_circ = QuantumCircuit(m,ne)
    
    for i in range(n)[::-1]:
        
        theta = 0
        
        for j in range(i + 1)[::-1]:
            k = i - j + 1
            
            if k > QFT_cutoff:
                break
            
            if b_bits[j]:
                theta += 2 * pi / (1<<k)
        
        if theta != 0:
            adder_circ.cp(theta, control_qubit, a_qubits[i])
            
    return adder_circ


def make_multiply_controlled_adder_circ(control_qubits, b_bits, a_qubits, m):
    adder_circ = QuantumCircuit(m,ne)
    
    for i in range(n)[::-1]:
        
        theta = 0
        
        for j in range(i + 1)[::-1]:
            k = i - j + 1
            
            if k > QFT_cutoff:
                break
            
            if b_bits[j]:
                theta += 2 * pi / (1<<k)
        
        if theta != 0:
            adder_circ.mcp(theta, control_qubits, a_qubits[i])
            
    return adder_circ



########################################




def modular_addition(a_qubits,ancillary_qubit,b,N,n,m):
    
    QFT_circ = make_QFT_circ(a_qubits,n,m)
    QFT_circ_inverse = QFT_circ.inverse()
    ones_complement_circ = ones_complement_circuit(a_qubits, m)
    
    # (a' + (N-b))' = a + b - N
    circuit = ones_complement_circ + QFT_circ + make_adder_circ(int_to_bin_LSB(N-b,n), a_qubits, m) + QFT_circ_inverse + ones_complement_circ
    
    # copy high / sign bit to the ancillary qubit
    circuit.cx(a_qubits[-1], ancillary_qubit)
    
    # a + b - N +? N = (a + b) % N
    circuit += QFT_circ + make_controlled_adder_circ(ancillary_qubit, int_to_bin_LSB(N,n), a_qubits, m) + QFT_circ_inverse
    
    # (a + b) % N - b
    circuit += ones_complement_circ + QFT_circ + make_adder_circ(int_to_bin_LSB(b,n), a_qubits, m) + QFT_circ_inverse + ones_complement_circ
    
    # uncompute the ancillary qubit
    circuit.x(a_qubits[-1])
    circuit.cx(a_qubits[-1], ancillary_qubit)
    circuit.x(a_qubits[-1])
    
    # (a + b) % N - b + b = (a + b) % N
    circuit += QFT_circ + make_adder_circ(int_to_bin_LSB(b,n), a_qubits, m) + QFT_circ_inverse
    
    return circuit

def controlled_modular_addition(control_qubits,a_qubits,ancillary_qubit,b,N,n,m):
    if type(control_qubits) is int:
        control_qubits = [control_qubits]
    
    QFT_circ = make_QFT_circ(a_qubits,n,m)
    QFT_circ_inverse = QFT_circ.inverse()
    ones_complement_circ = ones_complement_circuit(a_qubits, m)
    
    # (a' + (N-b))' = a + b - N
    circuit = ones_complement_circ + QFT_circ + make_multiply_controlled_adder_circ(control_qubits, int_to_bin_LSB(N-b,n), a_qubits, m) + QFT_circ_inverse + ones_complement_circ
    
    # copy high / sign bit to the ancillary qubit
    circuit.cx(a_qubits[-1], ancillary_qubit)
    
    # a + b - N +? N = (a + b) % N
    circuit += QFT_circ + make_multiply_controlled_adder_circ(control_qubits + [ancillary_qubit], int_to_bin_LSB(N,n), a_qubits, m) + QFT_circ_inverse
    
    # (a + b) % N - b
    circuit += ones_complement_circ + QFT_circ + make_multiply_controlled_adder_circ(control_qubits, int_to_bin_LSB(b,n), a_qubits, m) + QFT_circ_inverse + ones_complement_circ
    
    # uncompute the ancillary qubit
    circuit.x( a_qubits[-1])
    circuit.cx(a_qubits[-1], ancillary_qubit)
    circuit.x( a_qubits[-1])
    
    # (a + b) % N - b + b = (a + b) % N
    circuit += QFT_circ + make_multiply_controlled_adder_circ(control_qubits, int_to_bin_LSB(b,n), a_qubits, m) + QFT_circ_inverse
    
    circuit.mcx(control_qubits, ancillary_qubit)
    circuit.x(ancillary_qubit)
    
    

    return circuit


def modular_multiplication(x_qubits,y_qubits,ancillary_qubit,b,N,n,m):
    circuit = QuantumCircuit(m,ne)
    
    for i in range(n):
        circuit += controlled_modular_addition(x_qubits[i], y_qubits,ancillary_qubit,(b<<i)%N,N,n,m)
        
    return circuit

    
def conditional_in_place_modular_multiplication(control_qubit, x_qubits, y_qubits, ancillary_qubit ,b,N,n,m):
    # only support singly-controlled multiplication
    assert type(control_qubit) is int
        
    circuit = QuantumCircuit(m,ne)
    
    for i in range(n):
        circuit += controlled_modular_addition([x_qubits[i],control_qubit], y_qubits, ancillary_qubit, (b<<i)%N, N, n, m)
        
        
    circuit2 = QuantumCircuit(m,ne)
    b2 = mod_inverse(b, N)
    for i in range(n):
        circuit2 += controlled_modular_addition([y_qubits[i],control_qubit], x_qubits, ancillary_qubit, (b2<<i)%N, N, n, m)
        
    circuit += circuit2.inverse()
    
    # conditional swap.
    for i,j in zip(x_qubits, y_qubits):
        circuit.cswap(control_qubit, i, j)
    

    return circuit


"""
takes |z>|0> -> |z>|a^z mod N>

z stores exponent bits
x stores result
y is ancilla
"""
def modular_exponentiation(z_qubits, x_qubits, y_qubits, ancillary_qubit, a, N, n, m):
    circuit = QuantumCircuit(m,ne)
    
    a_tmp = a
    for z_qubit in z_qubits:
        circuit += conditional_in_place_modular_multiplication(z_qubit, x_qubits, y_qubits, ancillary_qubit, a_tmp, N, n, m)
        a_tmp = pow(a_tmp,2,N)
        
    return circuit
        


########################################

circuit =  QuantumCircuit(m,ne)

# initialize multiplication register to 1
x_bits = int_to_bin_LSB(1, n)
circuit += make_initialization_circuit(x_qubits, x_bits, m)


# hadamard the exponent bits
for z_qubit in z_qubits:
    circuit.h(z_qubit)


# modular exponentiation
circuit += modular_exponentiation(z_qubits, x_qubits, y_qubits, ancillary_qubit, a, N, n, m)


# finish with the inverse QFT  
# reverse bits b/c big vs little endian
circuit += make_QFT_circ(z_qubits[::-1],ne,m).inverse()

# and finally measure the output
# reverse bits b/c big vs little endian
circuit.measure(z_qubits[::-1], list(range(ne)))


########################################

print("[+] generated circuit")
with open("circuit.qasm","w") as f: f.write(circuit.qasm())
print("[+] saved circuit\n")


# can load circuit from str.
#qc = QuantumCircuit.from_qasm_str(qasm_str)


########################################


# simulate circuit

# Select the QasmSimulator from the Aer provider
simulator = Aer.get_backend('qasm_simulator')

# Execute and get counts
result = execute(circuit, simulator, shots=1024*8).result()
counts = result.get_counts(circuit)


########################################

# analyze results


# find order of a classically

powers_of_a = [1]
a_tmp = 1
for i in range(N):
    a_tmp = (a_tmp * a) % N
    powers_of_a.append(a_tmp)
    
order = 1 + powers_of_a[1:].index(1)



print("[+] order = {}".format(order))

t  = pow(a, order//2, N)
if (order % 2 == 0) and t != N-1:
    
    print("[+] recovered factors of N = {}".format(N))
    print("[+] {}\t{}".format(
        gcd(N, t - 1),
        gcd(N, t + 1) ))


########################################


# process the quantum results

q = 1 << ne
measurement_probabilities = np.zeros(q)
for i,j in counts.items():
    i = int(i,2)
    measurement_probabilities[i] = j
    
measurement_probabilities /= np.sum(measurement_probabilities)

plt.plot(measurement_probabilities)
plt.xlabel("QFT^-1 {exponent}")
plt.ylabel("probability")
plt.title(f"Shor's Algorithm Measurement Probabilities\nN={N}, a={a}, ne={ne}")
plt.savefig(f"output_distribution_N={N},a={a},ne={ne}.png", dpi=300)
plt.show()


print("")
for i in range(10):
    s = np.random.choice(np.arange(q), p=measurement_probabilities)
    f = Fraction(s/q).limit_denominator(N)
    recovered_order = f.denominator
    print("[+] s = {}\tf = {: <10}order = {}".format(s, str(f), recovered_order))



from math import gcd, log2, pi
import random
import progressbar

import multiprocessing
import sys

# Import Qiskit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import Aer, execute


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


def shuffle(x):
    return random.sample(x, len(x))

########################################

# circuit parameters

# number of bits + padding to store N
n = 64
# number of bits in exponent
ne = n * 2

# total number of qubits used
m = 2*n + ne + 1

# stores result
x_qubits = list(range(ne, n+ne))
# auxillary bits used for controlled multiplication
y_qubits = list(range(n+ne, 2*n + ne))
# exponent bits
z_qubits = list(range(ne))
# auxillary bit used in modular addition
ancillary_qubit = m - 1



num_pad_zeros = 2
N = 1149785239057719001

# base of exponent
a = 395318069365307

print("a = {}\nN = {}\n".format(a,N))


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
    
    for i in progressbar.progressbar(range(n)):
        circuit += controlled_modular_addition([x_qubits[i],control_qubit], y_qubits, ancillary_qubit, (b<<i)%N, N, n, m)
        
        
    circuit2 = QuantumCircuit(m,ne)
    b2 = mod_inverse(b, N)
    for i in progressbar.progressbar(range(n)):
        circuit2 += controlled_modular_addition([y_qubits[i],control_qubit], x_qubits, ancillary_qubit, (b2<<i)%N, N, n, m)
        
    circuit += circuit2.inverse()
    
    # conditional swap.
    for i,j in zip(x_qubits, y_qubits):
        circuit.cswap(control_qubit, i, j)
    

    return circuit


"""
qiskit wasn't designed to support circuits with billions of gates,
and the "circuit to qasm string" function is really slow
for large circuits (seems to take quadratic time for some reason).

so I only work with relatively short qiskit circuits at a time,
and concatenate the exported qasm together.
"""
def conditional_in_place_modular_multiplication_to_string(control_qubit, x_qubits, y_qubits, ancillary_qubit ,b,N,n,m):
    # only support singly-controlled multiplication
    assert type(control_qubit) is int
    
    res = ""
            
    for i in progressbar.progressbar(range(n)):
        circuit = controlled_modular_addition([x_qubits[i],control_qubit], y_qubits, ancillary_qubit, (b<<i)%N, N, n, m)
        s = circuit.qasm()
        s = s[find_nth(s, '\n', 4) + 1:]
        res += s
        
        
    b2 = mod_inverse(b, N)
    # note the reversed order. this is b/c we're doing the inverse circuit.
    for i in progressbar.progressbar(range(n)[::-1]):
        circuit = controlled_modular_addition([y_qubits[i],control_qubit], x_qubits, ancillary_qubit, (b2<<i)%N, N, n, m)
        
        s = circuit.inverse().qasm()
        s = s[find_nth(s, '\n', 4) + 1:]
        res += s
        
    
    # conditional swap.
    circuit = QuantumCircuit(m, ne)
    for i,j in zip(x_qubits, y_qubits):
        circuit.cswap(control_qubit, i, j)
        
    s = circuit.qasm()
    
    s = s[find_nth(s, '\n', 4) + 1:]
    res += s

    return res


def modular_exponentiation(e_qubits, x_qubits, y_qubits, ancillary_qubit, a_val, N, n, m):
    circuit = QuantumCircuit(m,ne)
    
    for e_qubit in e_qubits:
        circuit += conditional_in_place_modular_multiplication(e_qubit, x_qubits, y_qubits, ancillary_qubit, a_val, N, n, m)
        a_val = (a_val**2)%N
        
    return circuit


########################################

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start


circuit =  QuantumCircuit(m,ne)

# initialize multiplication register to 1
a_bits = int_to_bin_LSB(1, n)
circuit += make_initialization_circuit(x_qubits, a_bits, m)


# hadamard the exponent bits
for z_qubit in z_qubits:
    circuit.h(z_qubit)
    
    
base_directory = "circuit_parts"


with open(f"{base_directory}/circuit_initialization.qasm","w") as f: f.write(circuit.qasm())


# finish with the inverse QFT  
# reverse bits b/c big vs little endian
circuit = make_QFT_circ(z_qubits[::-1],ne,m).inverse()
circuit.measure(z_qubits[::-1], list(range(ne)))

with open(f"{base_directory}/circuit_finalization.qasm","w") as f:
    s = circuit.qasm()
    s = s[find_nth(s, '\n', 4):]
    f.write(s)


# filler text
s = "// not including the full circuit b/c it's absolutely massive\n"
for i in range(len(z_qubits)):
    with open(f"{base_directory}/circuit_z_{i}.qasm","w") as f:
        f.write(s)
    

# modular exponentiation
a_tmp = a
for i,z_qubit in enumerate(z_qubits[:]):
    s = conditional_in_place_modular_multiplication_to_string(z_qubit, x_qubits, y_qubits, ancillary_qubit, a_tmp, N, n, m)
    a_tmp = pow(a_tmp,2,N)
    
    with open(f"{base_directory}/circuit_z_{i}.qasm","w") as f:
        f.write(s)



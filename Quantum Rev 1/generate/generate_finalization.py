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


def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start



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



#QFT_circ = make_QFT_circ(a_qubits,n,m)
#QFT_circ_inverse = QFT_circ.inverse()
# QFT_circ.qasm()

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

#adder_circ = make_adder_circ(b_bits, a_qubits, m)

# multiple-control phase gate
# https://qiskit.org/textbook/ch-gates/more-circuit-identities.html#c-from-
# express in terms of rotations + toffoli, and decompose tofolli into 6 CNOTs


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



# circuit parameters

# 64-bit parameters
N = 183265954861256291
a = 74583193119493


# allocate registers
n = 64
m = ne = num_x_bits = 2*n

q_mem = [0 for _ in range(7*n + 2)]

z_qubits = [i for i in range(0, 2*n)]
b_qubits = [i for i in range(2*n, 3*n)]
N_qubits = [i for i in range(3*n, 4*n)]

tmp_qubits = [i for i in range(4*n, 5*n)]
c_qubits   = [i for i in range(5*n, 6*n)]
tmp2_qubits = [i for i in range(6*n, 7*n)]

ancilla_qubit = 7*n
high_qubit = 7*n + 1





########################################

    
    
base_directory = "circuit_parts"



# finish with the inverse QFT  
# reverse bits b/c big vs little endian
circuit = make_QFT_circ(z_qubits[::-1],ne,m).inverse()
circuit.measure(z_qubits[::-1], list(range(ne)))

with open(f"{base_directory}/circuit_finalization.qasm","w") as f:
    s = circuit.qasm()
    s = s[find_nth(s, '\n', 4) + 1:]
    f.write(s)







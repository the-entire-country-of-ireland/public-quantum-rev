import re
import progressbar 
from math import gcd, log2

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

"""
Curraco adder
https://arxiv.org/pdf/1905.09749.pdf
"""

def swap_gate(a,b):
    return f"swap q[{a}],q[{b}];\n"

def not_gate(a):
    return f"x q[{a}];\n"

def cnot_gate(a,b):
    return f"cx q[{a}],q[{b}];\n"

def ccnot_gate(a,b,c):
    return f"ccx q[{a}],q[{b}],q[{c}];\n"




# extract number within brackets
pattern = re.compile(r"\[(\d+)\]")

def eval_swap(i,j):
    q_mem[i],q_mem[j] = q_mem[j],q_mem[i]
    
def eval_ccx(i,j,k):
    if q_mem[i] and q_mem[j]:
        q_mem[k] = 1 - q_mem[k]
        
def eval_cx(i,j):
    if q_mem[i]:
        q_mem[j] = 1 - q_mem[j]
        
def eval_x(i):
    q_mem[i] = 1 - q_mem[i]

def eval_line(s):
    args = pattern.findall(s)
    args = list(map(int, args))
    
    if s.startswith("x "):
        eval_x(*args)
    elif s.startswith("cx "):
        eval_cx(*args)
    elif s.startswith("ccx "):
        eval_ccx(*args)
    elif s.startswith("swap "):
        eval_swap(*args)
        

"""
assume all the gates are self-inverse for now,
so we only have to reverse the order of the gates
"""
def invert_circuit(circuit):
    s = circuit.split("\n")
    s = s[:-1]
    s = s[::-1]
    s = '\n'.join(s)
    s += '\n'
    return s
        
#####################################

# we represent negative numbers via the ones-complement

# big-endian encoding (MSB first)
# ones complement
def int_to_bin(x,n):
    sign = (x < 0)
    x = abs(x)
    x = bin(x)[2:]
    x = x.zfill(n)[:n]
    if sign:
        x = x.replace("1","a").replace("0","1").replace("a","0")
    return x

def ones_complement(x,n):
    return (1<<n) + ~x


# MSB first
def bin_to_int(L):
    x = ''.join([str(i) for i in L])
    x = eval('0b' + x)
    return x

#####################################

# The in-place majority gate MAJ
def MAJ(a,b,c):
    s = ""
    s += cnot_gate(c,b)
    s += cnot_gate(c,a)
    s += ccnot_gate(a,b,c)
    
    return s

# UnMajority and Add
def UMA(a,b,c):
    s = ""
    s += ccnot_gate(a,b,c)
    s += cnot_gate(c,a)
    s += cnot_gate(a,b)
    
    return s

# addition circuit
"""
in-place addition of a + b
result is stored in `b`.
`a` is left unchanged.
a,b qubits must be specified in little endian order (LSB first)

`ancilla_qubit` must be initialized to 0; is reset to 0 at end of circuit

if given, then the high bit is stored in `high_qubit`
"""
def addition_circuit(n, a_qubits, b_qubits, ancilla_qubit, high_qubit=None):
    assert len(a_qubits) == n
    assert len(b_qubits) == n
    
    qubit_order = [ancilla_qubit]
    for i,j in zip(b_qubits, a_qubits):
        qubit_order += [i,j]
    
    s = ""
    for i in range(n):
        idxs = qubit_order[2*i:3+2*i]
        s += MAJ(*idxs)
    
    if high_qubit is not None:
        s += cnot_gate(a_qubits[-1], high_qubit)
        
    for i in range(n)[::-1]:
        idxs = qubit_order[2*i:3+2*i]
        s += UMA(*idxs)
        
    return s

def ones_complement_circuit(a_qubits):
    circuit = ""
    for i in a_qubits:
        circuit += not_gate(i)
    return circuit
    
# (a - b) = (a' + b)'
def subtraction_circuit(n, a_qubits, b_qubits, ancilla_qubit, high_qubit=None):
    
    circuit = ""
    circuit += ones_complement_circuit(a_qubits)
    
    circuit += addition_circuit(n, a_qubits, b_qubits, ancilla_qubit, high_qubit)
    
    circuit += ones_complement_circuit(a_qubits)
    circuit += ones_complement_circuit(b_qubits)
        
    return circuit

# (-a + b) = (a + b')'
def subtraction_circuit_2(n, a_qubits, b_qubits, ancilla_qubit, high_qubit=None):
    
    circuit = ""
    circuit += ones_complement_circuit(b_qubits)
    
    circuit += addition_circuit(n, a_qubits, b_qubits, ancilla_qubit, high_qubit)
    
    circuit += ones_complement_circuit(b_qubits)
        
    return circuit





"""
computes in-place modular addition:
|a>|b>|N> -> |a>|(a+b)%N>|N>

uses several ancilla qubits, which are all uncomputed to zero by the circuit.
"""
def modular_addition_circuit_(n, a_qubits, b_qubits, N_qubits, tmp_qubits, ancilla_qubit, high_qubit):
    
    # |a>|b>|N>|0>|0>
    
    circuit = addition_circuit(n, a_qubits, b_qubits, ancilla_qubit, None)
    # |a>|a+b>|N>|0>|0>
    
    circuit += subtraction_circuit_2(n, N_qubits, b_qubits, ancilla_qubit, high_qubit)
    # |a>|a+b-N>|N>|high qubit>|0>
    
    for i,j in zip(N_qubits, tmp_qubits):
        circuit += ccnot_gate(high_qubit, i, j)
    # |a>|a+b-N>|N>|high qubit>|0 or N>
    
    circuit += addition_circuit(n, tmp_qubits, b_qubits, ancilla_qubit, None)
    # |a>|(a+b)%N>|N>|high qubit>|0 or N>
    
    # uncompute tmp
    for i,j in zip(N_qubits, tmp_qubits):
        circuit += ccnot_gate(high_qubit, i, j)
    # |a>|(a+b)%N>|N>|high qubit>|0>
    
    # uncompute high qubit
    circuit += subtraction_circuit_2(n, a_qubits, b_qubits, ancilla_qubit, high_qubit)
    circuit += addition_circuit(n, a_qubits, b_qubits, ancilla_qubit, None)
    
    circuit += not_gate(high_qubit)
        
    return circuit

def memoize(func):
    cache = dict()
    
    def memoized_func(*args):
        x = []
        for i in args:
            if type(i) is list:
                i = tuple(i)
            x.append(i)
        x = tuple(x)
        args = x
    
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result

    return memoized_func

modular_addition_circuit = memoize(modular_addition_circuit_)

################################


# 256-bit parameters
# a = 3
# N = 177513599716362998539142178307007771385442861398395017247577542547366103

# 64-bit parameters
N = 183265954861256291
a = 74583193119493


# allocate registers
n = 64
num_x_bits = 2*n

q_mem = [0 for _ in range(7*n + 2)]

x_qubits = [i for i in range(0, 2*n)]
b_qubits = [i for i in range(2*n, 3*n)]
N_qubits = [i for i in range(3*n, 4*n)]

tmp_qubits = [i for i in range(4*n, 5*n)]
c_qubits   = [i for i in range(5*n, 6*n)]
tmp2_qubits = [i for i in range(6*n, 7*n)]

ancilla_qubit = 7*n
high_qubit = 7*n + 1

"""
if `conditional_qubit`, then set q_mem[dst] = src

Note: `src` contains the classical VALUES, not the index of any quantum registers.
"""
def conditional_classical_write(conditional_qubit, src_bits, dst_registers):
    s = ""
    for i,j in zip(src_bits, dst_registers):
        if i:
            s += cnot_gate(conditional_qubit, j)
            
    return s

def conditional_conditional_classical_write(conditional_qubit1, conditional_qubit2, src_bits, dst_registers):
    s = ""
    for i,j in zip(src_bits, dst_registers):
        if i:
            s += ccnot_gate(conditional_qubit1, conditional_qubit2, j)
            
    return s


"""
f_{a,N}(x) = (a*x)%N
for classical parameters a,N

a*x 
=
a * x_0 + 
(2*a) * x_1 +
(2^2*a) * x_2 +
...
"""
def conditional_multiplication(conditional_qubit, const_factors, b_qubits, a_qubits, c_qubits):
    circuit = ""
    
    for src_bits, b_bit in zip(const_factors, b_qubits):
    
        circuit += conditional_conditional_classical_write(conditional_qubit, b_bit, src_bits, a_qubits)
    
        circuit += modular_addition_circuit(n, a_qubits, c_qubits, N_qubits, tmp_qubits, ancilla_qubit, high_qubit)
        
        circuit += conditional_conditional_classical_write(conditional_qubit, b_bit, src_bits, a_qubits)
        
        
    # if the control is 0, then copy the argument into the dst register.
    circuit += not_gate(conditional_qubit)
    for i,j in zip(b_qubits, c_qubits):
        circuit += ccnot_gate(conditional_qubit, i, j)
    circuit += not_gate(conditional_qubit)
    
    return circuit



# multiplication

def multiplication(const_factors, b_qubits, a_qubits, c_qubits):
    circuit = ""
    for src_bits, b_bit in zip(const_factors, b_qubits):
    
        circuit += conditional_classical_write(b_bit, src_bits, a_qubits)
    
        circuit += modular_addition_circuit(n, a_qubits, c_qubits, N_qubits, tmp_qubits, ancilla_qubit, high_qubit)
        
        circuit += conditional_classical_write(b_bit, src_bits, a_qubits)

    return circuit


# classical constant factors, in little endian (LSB first)
def get_const_factors(a,N,n):
    const_factors = [(a<<i)%N for i in range(n)]
    const_factors = [int_to_bin(t,n)[::-1] for t in const_factors]
    const_factors = [[int(i) for i in t] for t in const_factors]

    return const_factors


def in_place_conditional_multiplication(a, N, num_x_bits, conditional_qubit, b_qubits, tmp2_qubits, c_qubits):
    circuit = ""
    circuit += conditional_multiplication(conditional_qubit, get_const_factors(a,N,n), b_qubits, tmp2_qubits, c_qubits)
    
    for i,j in zip(b_qubits, c_qubits):
        circuit += swap_gate(i,j)
    
    circuit += invert_circuit(conditional_multiplication(conditional_qubit, get_const_factors(mod_inverse(a, N),N,n), b_qubits, tmp2_qubits, c_qubits))
    
    return circuit



def in_place_modular_exponentiation(
        a,N,num_x_bits,
        x_qubits, b_qubits, tmp2_qubits, c_qubits):
    circuit = ""
    
    a_val = a
    for x_qubit in x_qubits:
        circuit += in_place_conditional_multiplication(a_val,N,num_x_bits,x_qubit, b_qubits, tmp2_qubits, c_qubits)
        a_val = (a_val**2)%N
        
    return circuit


################################



# create circuit
base_directory = "circuit_parts"


circuit = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[{}];
creg c[{}];

""".format(len(q_mem), len(x_qubits))

for i,v in enumerate(N_qubits):
    circuit += "cx c[{}],q[{}];\n".format(i,v)


# hadamard the exponent bits
for x_qubit in x_qubits:
    circuit += f"h q[{x_qubit}];\n"

# initialize `b` to 1
circuit += not_gate(b_qubits[0])


with open(f"{base_directory}/circuit_initialization.qasm","w") as f:
    f.write(circuit)




"""
Instead of writing the entire circuit to a file, we write each modular multiplication to its own circuit_part file.
This means that teams don't have to analyze the entire large circuit.
"""

# circuit += in_place_modular_exponentiation( a,N,num_x_bits, x_qubits, b_qubits, tmp2_qubits, c_qubits)


a_val = a
for i,x_qubit in enumerate(progressbar.progressbar(x_qubits)):
    circuit = in_place_conditional_multiplication(a_val,N,num_x_bits,x_qubit, b_qubits, tmp2_qubits, c_qubits)
    a_val = (a_val**2)%N
    
    with open(f"{base_directory}/circuit_part_{i}.qasm","w") as f:
        f.write(circuit)
        



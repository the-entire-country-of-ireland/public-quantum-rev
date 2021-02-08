import re
import math


# find nth occurence of needle in a haystack

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start


####################################################


with open("../challenge/circuit_parts/circuit_part_0.qasm", "r") as f:
    prog = f.read()
    

####################################################

# extract the 1st addition subroutine

s = "cx q[384],q[320];\n"
i = prog.index(s)
addition_circuit = prog[i:]

s = "cx q[448],q[320];\n"
j = addition_circuit.index(s) + len(s)
addition_circuit = addition_circuit[:j]

len_addition_circuit = addition_circuit.count("\n")

"""
note that all the additions and registers and whatnot are 64-bit integers.
and 384 / 64 = 6, 320 / 64  = 5,
so this is adding-in-place registers 5 and 6
"""

prog = prog.replace(addition_circuit, "add #5, #6\n")

####################################################

# coalesce all the repeated x-gates

s = "x q[320];\n"
i = prog.index(s)
negation_circuit = prog[i:]

s = "x q[383];\n"
j = negation_circuit.index(s) + len(s)
negation_circuit = negation_circuit[:j]

prog = prog.replace(negation_circuit, "negate #5\n")

####################################################

"""
once you've reversed the first addition subroutine,
then identifying the rest of them should be pretty straightforward.
"""

# extract the 2nd addition subroutine

s = "cx q[192],q[320];\n"
i = prog.index(s)
addition_circuit = prog[i:]

j = find_nth(addition_circuit, "\n", len_addition_circuit + 1) + 1
addition_circuit = addition_circuit[:j]

#this is adding-in-place registers 5 and 3
prog = prog.replace(addition_circuit, "add #5, #3\n")

####################################################

# coalesce the ccx

s = "ccx q[449],q[192],q[64];\n"
i = prog.index(s)
negation_circuit = prog[i:]

s = "ccx q[449],q[255],q[319];\n"
j = negation_circuit.index(s) + len(s)
negation_circuit = negation_circuit[:j]

prog = prog.replace(negation_circuit, "ccx q[449], #3, #4\n")

####################################################

# extract the 3rd addition subroutine

s = "cx q[64],q[320];\n"
i = prog.index(s)
addition_circuit = prog[i:]

j = find_nth(addition_circuit, "\n", len_addition_circuit + 1) + 1
addition_circuit = addition_circuit[:j]

prog = prog.replace(addition_circuit, "add #4, #5\n")

####################################################

# extract the 4th addition subroutine

s = "cx q[384],q[320];\n"
i = prog.index(s)
addition_circuit = prog[i:]

j = find_nth(addition_circuit, "\n", len_addition_circuit + 1) + 1
addition_circuit = addition_circuit[:j]

prog = prog.replace(addition_circuit, "add #6, #5\n")


####################################################

"""
coalesce ccx writes
"""

N = 64
# register #6 starts with this qubit
idx_start = 6 * N


def bin_to_int(b):
    res = 0
    for i in b:
        res |= 1<<i
    return res

"""
match repeated ccx lines where the first two controls are the same

ccx q[0],q[128],q[384];
ccx q[0],q[128],q[386];
...
"""
pattern = re.compile(r"ccx q\[(\d+)\],q\[(\d+)],q\[\d+\];\n(ccx q\[\1\],q\[\2],q\[\d+\];\n)+")

pattern2 = re.compile(r"ccx q\[(\d+)\],q\[(\d+)],q\[(\d+)\];")

indices = [(m.start(0), m.end(0)) for m in re.finditer(pattern, prog)]


new_prog = prog

for segment in indices:
    s = prog[segment[0]:segment[1]]
    
    qubits = pattern2.findall(s)
    qubits = [int(i[2]) for i in qubits]
    
    idx_min = min(qubits)
    idx_max = max(qubits)
    register = math.floor(idx_min / N)
    
    if idx_max >= register * N + N:
        continue
    
    qubits = [i - register*N for i in qubits]
    
    number = bin_to_int(qubits)
    
    i = find_nth(s, ",q", 2)
    cmd = s[:i]
    cmd += f", `{number}` -> #{register}\n"
    
    new_prog = new_prog.replace(s, cmd)
    
prog = new_prog


####################################################


print(prog[:prog.index("swap")])


import re
from math import pi, log2
import numpy as np

# grep mcphase ../challenge/circuit_parts/circuit_z_0.qasm -A 0 > output_z_0.txt

with open("output_z_0.txt", "r")  as f:
    data = f.read(10000000)
    
data = data.split("\n--\n")


n = 64

pattern = re.compile(r"\[(\d+)\]")

print("\n########################################\n")


lines = data[0].split("\n")

phases = np.zeros(n)

for line in lines:
    args = pattern.findall(line)
    i = int(args[-1]) - 3 * n
    
    phase = line[line.find("(")+1:line.find(")")]
    phase = eval(phase) # b/c we love security and best practices here
    print(phase / (2 * pi))
    
    phases[i] = phase / (2 * pi)
    

print("\n########################################\n")


"""
we need to convert the rotation angles into binary decimals.
Let's sweep the number of digits and see what the resulting loss of accuracy is.
"""
for m in range(5, 20):
    scaling = (1 << m)
    scaled_phases = phases * scaling
    error = np.sum(np.abs(np.round(scaled_phases) - scaled_phases))
    print("[+] {}\t{}".format(m, error))
    

print("\n########################################\n")


m = 12
scaling = (1 << m)
scaled_phases = [int(i) for i in np.round(phases * scaling)]

bits = []
for i,v in enumerate(scaled_phases):
    print(f"[+] {i}\t{v:0{m}b}")
    bits.append(f"{v:0{m}b}"[0])
    
number = int(''.join(bits), 2)


print("\n########################################\n")


def process(lines, m=12):
    phases = np.zeros(n)
    
    for line in lines:
        args = pattern.findall(line)
        i = int(args[-1]) - 3 * n
        
        phase = line[line.find("(")+1:line.find(")")]
        phase = eval(phase) # b/c we love security and best practices here
        
        phases[i] = phase / (2 * pi)
    
    scaling = (1 << m)
    scaled_phases = [int(i) for i in np.round(phases * scaling)]
    
    bits = []
    for i,v in enumerate(scaled_phases):
        bits.append(f"{v:0{m}b}"[0])
    
    number = int(''.join(bits[::-1]), 2)
    return number


for d in data[:20]:
    lines = d.split("\n")
    number = process(lines)
    print(number)


N = process(data[1].split("\n"))
a = process(data[2].split("\n"))

print(f"\n\nrecovered N = {N}")

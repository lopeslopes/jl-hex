import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


angles = []
AA = []
BA = []
AB = []
BB = []
with open("AAstack_separations") as f:
    data = f.readlines()
    for line in data:
        if line[0:2] == "0.":
            angles.append(float(line))
        elif line[0:3] == "AA:":
            aux = line.split(":")
            AA.append(float(aux[1]))
        elif line[0:3] == "BA:":
            aux = line.split(":")
            BA.append(float(aux[1]))
        elif line[0:3] == "AB:":
            aux = line.split(":")
            AB.append(float(aux[1]))
        elif line[0:3] == "BB:":
            aux = line.split(":")
            BB.append(float(aux[1]))


ax1 = plt.subplot()
ax1.plot(angles, AA)
ax1.plot(angles, BA)
ax1.plot(angles, AB)
ax1.plot(angles, BB)

plt.legend(["AA", "BA", "AB", "BB"])
plt.ylim((0,0.03))


plt.show()

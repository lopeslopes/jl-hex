import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


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


min_index = np.argmin(AB)
print("Angle:      ", angles[min_index])
print("Separation: ", AB[min_index])


ax1 = plt.subplot()
ax1.plot(angles, AA, color="c")
# ax1.plot(angles, BA, color="orange")
ax1.plot(angles, AB, color="green")
# ax1.plot(angles, BB, color="red")

for q in [61.0, 60.0, 59.0, 58.0]:
    p = 1.0
    magic_angle = np.acos((3.0*(q**2) - (p**2))/(3.0*(q**2) + (p**2)))
    ax1.axvline(magic_angle, color="black")

plt.xlabel(r"Angle (rad)")
plt.ylabel(r"Separation ($\AA$)")

plt.legend(["AA", "AB"])
# plt.legend(["AA", "BA", "AB", "BB"])
plt.ylim((0,0.03))


plt.show()

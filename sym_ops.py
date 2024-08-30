import numpy as np
import spglib
import matplotlib.pyplot as plt
from sympy import KroneckerDelta as delta
from sympy import LeviCivita as lciv


def read_lattice(path):
    try:
        with open(file=path) as f:
            data = f.readlines()
        points = []
        for line in data:
            aux = line.split(";")
            points.append([float(aux[0]), float(aux[1]), 0.0])

        # points = np.array(points)
        return points
    except:
        print("No file found in "+path)
        return 0


a = 2.46
d = np.sqrt((a**2) / (2 * (1 - np.cos(2 * np.pi / 3))))
d1 = [d * np.cos(np.pi / 6), d * np.sin(np.pi / 6), 0.0]

# DEFINE LATTICE VECTORS

# HEXAGONAL CELL MATCHING OVERLAP LATTICES
lattice = np.array([[4428.000519758899, 2484.072918780167, 0.0],
                    [-4365.2705123956675, 2592.7244786925335, 0.0],
                    [0.0, 0.0, 10000.0]])

# SIMPLE CUBIC LATTICE FOR TESTING
# lattice = np.array([[1.0, 0.0, 0.0],
#                     [0.0, 1.0, 0.0],
#                     [0.0, 0.0, 1.0]])

# BCC LATTICE FOR TESTING
# lattice = np.array([[-0.5, 0.5, 0.5],
#                     [0.5, -0.5, 0.5],
#                     [0.5, 0.5, -0.5]])

# FCC LATTICE FOR TESTING
# lattice = np.array([[0.5, 0.5, 0.0],
#                     [0.0, 0.5, 0.5],
#                     [0.5, 0.0, 0.5]])

# TESTING LENGTHS
dist1 = np.sqrt(lattice[0,0]**2 + lattice[0,1]**2 + lattice[0,2]**2)
dist2 = np.sqrt(lattice[1,0]**2 + lattice[1,1]**2 + lattice[1,2]**2)
dist3 = np.sqrt(lattice[2,0]**2 + lattice[2,1]**2 + lattice[2,2]**2)
print(f"|a1| = {dist1}")
print(f"|a2| = {dist2}")
print(f"|a3| = {dist3}")

# TESTING OTHER ATOM POSITIONS TO DEBUG SPGLIB
# positions = np.array([[0.0, 0.0, 0.0],
#                       [-1322.249830067740, 2493.3053571691007, 0.0]])
# numbers = [1, 2]
#
# positions = np.array([[0.0, 0.0, 0.0],
#                       [-1322.249830067740, 2493.3053571691007, 0.0],
#                       [1496.9101757073777, 2520.290105417622, 0.0],
#                       [-1.2306879586268324, 128.5354838445569, 0.0]])
# numbers = [1, 2, 3, 4]

positions = np.array([[0.0, 0.0, 0.0]])
numbers = [1]

graphene_structure = (lattice, positions, numbers)

# SYMMETRY CALCULATIONS
sym_data = spglib.get_symmetry_dataset(graphene_structure, symprec=1e-5)
print(f"Group number:           {sym_data.number}")
print(f"Hall number:            {sym_data.hall_number}")
print(f"International notation: {sym_data.international}")
print(f"Hall notation:          {sym_data.hall}")

# READING LATTICES GENERATED BY JULIA PROGRAM
atomsAA = np.array(read_lattice("data/0.0191435_0.001/latticeAA.dat"))
lenAA = len(atomsAA)
atomsAB = np.array(read_lattice("data/0.0191435_0.001/latticeAB.dat"))
lenAB = len(atomsAB)
atomsBA = np.array(read_lattice("data/0.0191435_0.001/latticeBA.dat"))
lenBA = len(atomsBA)
atomsBB = np.array(read_lattice("data/0.0191435_0.001/latticeBB.dat"))
lenBB = len(atomsBB)

# PLOTTING REAL SPACE POINTS AFTER ROTATIONS
ax2 = plt.subplot(111, projection="3d")

p1 = lattice[0,:]
p2 = lattice[1,:]
p3 = lattice[2,:]

print(f"Number of symmetry operations of group: {len(sym_data.rotations)}")

for j in range(1):
    new_points = []
    spin = []
    grp_chr = []
    op_viz = []
    pn = lattice[j,:]
    for i, rot in enumerate(sym_data.rotations):
        if (1):
            # DISSECATING ROTATION MATRICES
            trc = np.trace(rot)
            values, vectors = np.linalg.eig(rot)
            det = np.linalg.det(rot)
            trc_det = trc * det

            if (trc_det == 3):
                gen_rot = rot

            else:
                if (trc_det == -1):
                    # IMPROPER ROTATION
                    if (trc == 1):
                        ax_ind = np.where(values == -1.0)[0][0]
                        rot_axis = vectors[ax_ind]
                        rot_angle = 0.0

                    # PROPER ROTATION
                    elif (trc == -1):
                        ax_ind = np.where(values == 1.0)[0][0]
                        rot_axis = vectors[ax_ind]
                        rot_angle = np.pi

                else:
                    if (det == 1.0):
                        ax_ind = np.where(values == 1.0)[0][0]
                        rot_axis = vectors[ax_ind]
                        rot_angle = np.acos((trc - 1)/2)

                    elif (det == -1.0):
                        ax_ind = np.where(values == -1.0)[0][0]
                        rot_axis = vectors[ax_ind]
                        rot_angle = np.acos((trc + 1)/2)

                gen_rot = np.zeros((3,3), dtype=float)
                max_viz = 15
                for k in range(max_viz+1):
                    alpha_rot = rot_angle*float(k/max_viz)
                    for m in range(3):
                        for l in range(3):
                            aux_set = set([0, 1, 2])
                            aux_set.discard(m)
                            aux_set.discard(l)
                            n = aux_set.pop()
                            gen_rot[m,l] = delta(m,l)*np.cos(alpha_rot) + (det-np.cos(alpha_rot))*rot_axis[m]*rot_axis[l] - np.sin(alpha_rot)*lciv(m,l,n)*rot_axis[n]

                    dgd_point = gen_rot @ pn
                    op_viz.append(dgd_point.real)

                op_viz = np.transpose(np.array(op_viz))
                ax2.plot(op_viz[0,:], op_viz[1,:], op_viz[2,:])
                op_viz = []

            # POINT
            new_p = gen_rot @ pn
            new_points.append(new_p.real)

            # # SPIN TEST
            # spin_test = pn + [0.0, 0.0, 1.0]
            # new_spin_test = gen_rot.real @ spin_test
            # ds = new_spin_test - new_p
            # spin.append(int(ds[2].real))
            # gc = 0
            # if all(np.isclose(new_p, pn)):
            #     gc += int(ds[2])
            # grp_chr.append(gc)

    new_points = np.transpose(np.array(new_points))
    ax2.scatter(new_points[0,:], new_points[1,:], new_points[2,:])

# print(grp_chr)

ax2.quiver(0.0, 0.0, 0.0, lattice[0,0], lattice[0,1], 0.0, color="red")
ax2.quiver(0.0, 0.0, 0.0, lattice[1,0], lattice[1,1], 0.0, color="green")

ax2.set_box_aspect((1, 1, 1))
ax2.set_xlim([-10000.0, 10000.0])
ax2.set_ylim([-10000.0, 10000.0])
ax2.set_zlim([-10000.0, 10000.0])

plt.show()

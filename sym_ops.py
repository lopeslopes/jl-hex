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
                    [0.0, 0.0, 10000.0]], dtype=np.float64)

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

positions = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
numbers = [1]

graphene_structure = (lattice, positions, numbers)

# SYMMETRY CALCULATIONS
sym_data = spglib.get_symmetry_dataset(graphene_structure, symprec=1e-5)
print(f"Group number:           {sym_data.number}")
print(f"Hall number:            {sym_data.hall_number}")
print(f"International notation: {sym_data.international}")
print(f"Hall notation:          {sym_data.hall}")

# READING LATTICES GENERATED BY JULIA PROGRAM
atomsAA = np.array(read_lattice("data/0.0191435_0.001/latticeAA.dat"), dtype=np.float64)
lenAA = len(atomsAA)
atomsAB = np.array(read_lattice("data/0.0191435_0.001/latticeAB.dat"), dtype=np.float64)
lenAB = len(atomsAB)
atomsBA = np.array(read_lattice("data/0.0191435_0.001/latticeBA.dat"), dtype=np.float64)
lenBA = len(atomsBA)
atomsBB = np.array(read_lattice("data/0.0191435_0.001/latticeBB.dat"), dtype=np.float64)
lenBB = len(atomsBB)

# PLOTTING REAL SPACE POINTS AFTER ROTATIONS
ax2 = plt.subplot(111, projection="3d")

p1 = lattice[0,:]
p2 = lattice[1,:]
p3 = lattice[2,:]

print(f"Number of symmetry operations of group: {len(sym_data.rotations)}")

# MAKING LIST OF POINTS IN HEXAGON VERTICES
def_positions = np.transpose([lattice[0,:],
                              lattice[1,:]*dist1/dist2,
                              -lattice[0,:],
                              -lattice[1,:]*dist1/dist2,
                              [0.0, dist1, 0.0],
                              [0.0, -dist1, 0.0]])

ax2.scatter(def_positions[0,:], def_positions[1,:], def_positions[2,:], s=200)
def_positions = np.transpose(def_positions)

aux_vec1 = lattice[0,:]/np.linalg.norm(lattice[0,:])
aux_vec2 = lattice[1,:]/np.linalg.norm(lattice[1,:])
aux_vec3 = np.array([0.0, 1.0, 0.0]) + aux_vec1
aux_vec3 = aux_vec3/np.linalg.norm(aux_vec3)
aux_vec4 = np.array([0.0, 1.0, 0.0]) + aux_vec2
aux_vec4 = aux_vec4/np.linalg.norm(aux_vec4)

grp_chr = np.zeros(len(sym_data.rotations), dtype=int)
grp_chr_names = ["" for i in range(len(sym_data.rotations))]
for i, rot in enumerate(sym_data.rotations):
    # DISSECATING ROTATION MATRICES
    values, vectors = np.linalg.eig(rot)
    vectors = np.transpose(vectors)
    trc = np.trace(rot)
    det = np.linalg.det(rot)
    trc_det = trc * det

    if (det == 1.0):
        if (trc == -1):
            ax_ind = np.where(values == 1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            if all(np.isclose(rot_axis, np.array([0.0, 0.0, 1.0]))):
                grp_chr_names[i] = "C2"
            elif all(np.isclose(rot_axis, np.array([1.0, 0.0, 0.0]))):
                grp_chr_names[i] = "C*2"
            elif all(np.isclose(rot_axis, np.array([0.0, 1.0, 0.0]))):
                grp_chr_names[i] = "C**2"
            elif (rot_axis[0] == rot_axis[1]):
                rot_axis = aux_vec3
                grp_chr_names[i] = "C*2"
            elif (rot_axis[0] == -rot_axis[1]):
                rot_axis = aux_vec4
                grp_chr_names[i] = "C*2"
            elif (rot_axis[0] == 2*rot_axis[1]):
                rot_axis = aux_vec1
                grp_chr_names[i] = "C**2"
            elif (2*rot_axis[0] == rot_axis[1]):
                rot_axis = aux_vec2
                grp_chr_names[i] = "C**2"
        elif (trc == 0):
            ax_ind = np.where(values == 1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            aux_name = "C3"
            if aux_name in grp_chr_names:
                rot_angle = rot_angle + 2*(np.pi - rot_angle)
            grp_chr_names[i] = "C3"
        elif (trc == 2):
            ax_ind = np.where(values == 1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            aux_name = "C6"
            if aux_name in grp_chr_names:
                rot_angle = rot_angle + 2*(np.pi - rot_angle)
            grp_chr_names[i] = "C6"
        elif (trc == 3):
            ax_ind = np.where(values == 1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            grp_chr_names[i] = "E"
        else:
            rot_axis = [0.0, 0.0, 1.0]
            rot_angle = 0.0

    elif (det == -1.0):
        if (trc == 1):
            ax_ind = np.where(values == -1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            if all(np.isclose(rot_axis, np.array([0.0, 0.0, 1.0]))):
                grp_chr_names[i] = "sigma_h"
            elif all(np.isclose(rot_axis, np.array([1.0, 0.0, 0.0]))):
                grp_chr_names[i] = "sigma_v"
            elif all(np.isclose(rot_axis, np.array([0.0, 1.0, 0.0]))):
                grp_chr_names[i] = "sigma_d"
            elif (rot_axis[0] == rot_axis[1]):
                rot_axis = aux_vec3
                grp_chr_names[i] = "sigma_v"
            elif (rot_axis[0] == -rot_axis[1]):
                rot_axis = aux_vec4
                grp_chr_names[i] = "sigma_v"
            elif (rot_axis[0] == 2*rot_axis[1]):
                rot_axis = aux_vec1
                grp_chr_names[i] = "sigma_d"
            elif (2*rot_axis[0] == rot_axis[1]):
                rot_axis = aux_vec2
                grp_chr_names[i] = "sigma_d"
        elif (trc == 0):
            ax_ind = np.where(values == -1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            aux_name = "S6"
            if aux_name in grp_chr_names:
                rot_angle = rot_angle + 2*(np.pi - rot_angle)
            grp_chr_names[i] = "S6"
        elif (trc == -2):
            ax_ind = np.where(values == -1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            aux_name = "S3"
            if aux_name in grp_chr_names:
                rot_angle = rot_angle + 2*(np.pi - rot_angle)
            grp_chr_names[i] = "S3"
        elif (trc == -3):
            ax_ind = np.where(values == -1.0)[0][0]
            rot_axis = vectors[ax_ind].real/np.linalg.norm(vectors[ax_ind].real)
            rot_angle = np.arccos((trc - det)/2)
            grp_chr_names[i] = "i"
        else:
            rot_axis = [0.0, 0.0, 1.0]
            rot_angle = 0.0

    gen_rot = np.zeros((3,3), dtype=np.float64)
    for m in range(3):
        for l in range(3):
            aux_set = set([0, 1, 2])
            aux_set.discard(m)
            aux_set.discard(l)
            n = aux_set.pop()
            gen_rot[m,l] = delta(m,l)*np.cos(rot_angle) + (det-np.cos(rot_angle))*rot_axis[m]*rot_axis[l] - np.sin(rot_angle)*lciv(m,l,n)*rot_axis[n]

    for jj, pn in enumerate(def_positions):
        new_points = []
        op_viz = []

        # # VISUALIZATION
        # max_viz = 100
        # for k in range(max_viz+1):
        #     alpha_rot = rot_angle*float(k/max_viz)
        #     for m in range(3):
        #         for l in range(3):
        #             aux_set = set([0, 1, 2])
        #             aux_set.discard(m)
        #             aux_set.discard(l)
        #             n = aux_set.pop()
        #             gen_rot[m,l] = delta(m,l)*np.cos(alpha_rot) + (det-np.cos(alpha_rot))*rot_axis[m]*rot_axis[l] - np.sin(alpha_rot)*lciv(m,l,n)*rot_axis[n]
        #
        #     dgd_point = gen_rot @ pn
        #     op_viz.append(dgd_point.real)
        #
        # op_viz = np.transpose(np.array(op_viz))
        # if (jj == 0):
        #     ax2.plot(op_viz[0,:], op_viz[1,:], op_viz[2,:])
        #
        # op_viz = []

        # POINT
        new_p = gen_rot @ pn
        new_points.append(new_p.real)

        # SPIN TESTING
        spin_test = pn + [0.0, 0.0, 1.0]
        new_spin_test = gen_rot.real @ spin_test
        ds = new_spin_test - new_p
        gc = 0
        if (new_p in def_positions):
            gc += int(ds[2])
        else:
            new_points.append(new_p.real)
        grp_chr[i] += gc

    new_points = np.transpose(np.array(new_points))
    ax2.scatter(new_points[0,:], new_points[1,:], new_points[2,:])

names_unique, counts = np.unique(grp_chr_names, return_counts=True)
chr_unique = np.zeros(len(names_unique), dtype=int)

for j, nu in enumerate(names_unique):
    ind = [i for i, val in enumerate(grp_chr_names) if val == nu]
    chr_final = 0
    cont = 0
    for i in ind:
        if grp_chr[i] != chr_final:
            # print(nu, chr_final, grp_chr[i])
            chr_final = grp_chr[i]
            cont += 1
    chr_unique[j] = chr_final

print("----------------------------")
print("CHARACTERES")
print(f"{"sym_op": ^7}", f"{"chr": ^3}", f"{"mult": ^4}")
for n, c in zip(names_unique, chr_unique):
    mult = counts[np.where(names_unique == n)[0][0]]
    print(f"{n: ^7}", f"{c: ^3}", f"{mult: ^4}")

# ax2.quiver(0.0, 0.0, 0.0, lattice[0,0], lattice[0,1], 0.0, color="red")
# ax2.quiver(0.0, 0.0, 0.0, lattice[0,0], lattice[0,1], 0.0, color="red")

ax2.set_box_aspect((1, 1, 1))
ax2.set_xlim([-10000.0, 10000.0])
ax2.set_ylim([-10000.0, 10000.0])
ax2.set_zlim([-10000.0, 10000.0])

plt.show()

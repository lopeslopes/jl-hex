import numpy as np
import spglib
import matplotlib.pyplot as plt


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
# lattice = np.array([[a, 0.0, 0.0],
#                     [a/2, a*np.sqrt(3)/2, 0.0],
#                     [0.0, 0.0, 2000.0]])

# lattice = np.array([[2931.090344051522, -36.2171866374553, 0.0],
#                     [-1434.1801683441458, 2556.5072920550774, 0.0],
#                     [0.0, 0.0, 50000.0]])

lattice = np.array([[4428.000519758899, 2484.072918780167, 0.0],
                    [62.730007363230285, 5076.797397472701, 0.0],
                    [0.0, 0.0, 50000.0]])

# DEFINE ATOMS IN UNIT CELL

# atoms = [[0.0, 0.0, 0.0],
#         [a/2, a*np.sqrt(3)/2, 0.0],
#         [a, 0.0, 0.0],
#         d1]

atomsAA = read_lattice("data/0.0191435_0.001/latticeAA.dat")
lenAA = len(atomsAA)
numbersAA = np.ones(lenAA)

atomsAB = read_lattice("data/0.0191435_0.001/latticeAB.dat")
lenAB = len(atomsAB)
numbersAB = 2 * np.ones(lenAB)

atomsBA = read_lattice("data/0.0191435_0.001/latticeBA.dat")
lenBA = len(atomsBA)
numbersBA = 3 * np.ones(lenBA)

atomsBB = read_lattice("data/0.0191435_0.001/latticeBB.dat")
lenBB = len(atomsBB)
numbersBB = 4 * np.ones(lenBB)

positions = np.concatenate((atomsAA, atomsAB))
numbers = np.concatenate((numbersAA, numbersAB))

graphene_structure = (lattice, positions, numbers)

# SYMMETRY CALCULATIONS
sym_data = spglib.get_symmetry_dataset(graphene_structure, symprec=1e-5)
print(f"{sym_data.number=}")
print(f"{sym_data.hall_number=}")
print(f"{sym_data.international=}")
print(f"{sym_data.hall=}")

# GAMMA CENTERED MESH
mesh = [10, 10, 10]
mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, graphene_structure, is_shift=[0, 0, 0])
# for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
#    print("%3d ->%3d %s" % (i, ir_gp_id, gp.astype(float) / mesh))

# IRREDUCIBLE K-POINTS
print("Number of ir-kpoints: %d" % len(np.unique(mapping)))
unique_kpts = np.transpose(grid[np.unique(mapping)] / np.array(mesh, dtype=float))

atomsAA = np.array(atomsAA)
atomsBA = np.array(atomsBA)
atomsAB = np.array(atomsAB)
atomsBB = np.array(atomsBB)
# print(atomsAA)

# PLOT OF K POINTS
ax1 = plt.subplot(211, projection="3d")
ax1.scatter(unique_kpts[0, :], unique_kpts[1, :], unique_kpts[2, :])

# PLOT OF DIRECT SPACE THINGS
ax2 = plt.subplot(212)
ax2.scatter(atomsAA[:,0], atomsAA[:,1], s=20, color="c")
ax2.scatter(atomsBA[:,0], atomsBA[:,1], s=20, color="green")
ax2.scatter(atomsAB[:,0], atomsAB[:,1], s=20, color="orange")
ax2.scatter(atomsBB[:,0], atomsBB[:,1], s=20, color="red")

ax2.quiver(0.0, 0.0, lattice[0,0], lattice[0,1], angles="xy", scale_units="xy", scale=1)
ax2.quiver(0.0, 0.0, lattice[1,0], lattice[1,1], angles="xy", scale_units="xy", scale=1)

ax1.set_aspect("equal")
ax2.set_aspect("equal")

plt.show()

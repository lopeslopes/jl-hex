import numpy as np
import spglib
import matplotlib.pyplot as plt

a = 2.46
d = np.sqrt((a**2) / (2 * (1 - np.cos(2 * np.pi / 3))))
d1 = [d * np.cos(np.pi / 6), d * np.sin(np.pi / 6), 0.0]
lattice = np.array([[a, 0.0, 0.0],
                    [a/2, a*np.sqrt(3)/2, 0.0],
                    [0.0, 0.0, 2000.0]])  # a large c to simulate a 2D material

# Define basis atoms in the unit cell
atoms = [[0.0, 0.0, 0.0],
         [a/2, a*np.sqrt(3)/2, 0.0],
         [a, 0.0, 0.0],
         d1]
positions = np.array(atoms)
numbers = [1, 1, 1, 2]

graphene_structure = (lattice, positions, numbers)

print(spglib.get_symmetry_dataset(graphene_structure, symprec=1e-5))
mesh = [20, 20, 20]

# Gamma centre mesh
mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, graphene_structure, is_shift=[0, 0, 0])
#for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
#    print("%3d ->%3d %s" % (i, ir_gp_id, gp.astype(float) / mesh))

# Irreducible k-points
print("Number of ir-kpoints: %d" % len(np.unique(mapping)))
unique_kpts = np.transpose(grid[np.unique(mapping)] / np.array(mesh, dtype=float))
all_kpts = np.transpose(grid)

ax1 = plt.subplot(111, projection="3d")
ax1.scatter(unique_kpts[0,:], unique_kpts[1,:], unique_kpts[2,:])
#ax1.scatter(all_kpts[0,:], all_kpts[1,:], all_kpts[2,:])
plt.show()

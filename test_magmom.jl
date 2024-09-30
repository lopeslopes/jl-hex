include("hex_utils.jl")
using .HexUtils
using PyCall
using LinearAlgebra
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Combinatorics
using StatsBase

spglib = pyimport("spglib")
ase = pyimport("ase")
sg = pyimport("ase.spacegroup")
viz = pyimport("ase.visualize.plot")

lattice = [[4428.000519758899, 2484.072918780167, 0.0],
           [-4365.2705123956675, 2592.7244786925335, 0.0],
           [0.0, 0.0, 10000.0]]

# lattice = [ [sqrt(3)/2, 1/2, 0.0],
#             [-sqrt(3)/2, 1/2, 0.0],
#             [0.0, 0.0, 10000]]

a1 = lattice[1,:][1]
a2 = lattice[2,:][1]
a3 = lattice[3,:][1]

cell_vol = dot(a1, cross(a2, a3))
b1 = 2.0*pi*cross(a2, a3)/cell_vol
b2 = 2.0*pi*cross(a3, a1)/cell_vol
b3 = 2.0*pi*cross(a1, a2)/cell_vol

positions = [[0.0, 0.0, 0.0]]
numbers = [1]
graphene = (lattice, positions, numbers)
sym_data = spglib.get_symmetry_dataset(graphene)

mesh = [20, 20, 20]
mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, graphene, is_shift=[0, 0, 0])

un_kpoints = []
un_ind = unique(mapping)
for j in un_ind
    i = j+1
    kp = grid[i,1]*b1 .+ grid[i,2]*b2 .+ grid[i,3]*b3
    push!(un_kpoints, kp)
end
un_kpoints = hcat(un_kpoints...)
un_kpoints = transpose(un_kpoints)

ax1 = subplot(211, projection="3d")
ax1.scatter(un_kpoints[:,1], un_kpoints[:,2], un_kpoints[:,3], color="c")

# ASE FEATURES TESTING
supercell = ase.cell.Cell.new(transpose(lattice))
println(supercell.lengths())

kcell = supercell.reciprocal()
println(kcell.lengths())

# ax2 = subplot(212)
# graphene_test = sg.crystal("C", [(0,0,0)], spacegroup=183, cellpar=supercell.cellpar(), size=(5,5,5))
# viz.plot_atoms(graphene_test, ax2)

# path = supercell.bandpath("KGLK", 100)
# test_path = path.cartesian_kpts()
# ax1.scatter(test_path[:,1], test_path[:,2], test_path[:,3], color="red")

ax1.set_aspect("equal")
show()

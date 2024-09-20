include("hex_utils.jl")
using .HexUtils
using PyCall
using PyCall: LinearAlgebra
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Combinatorics
using StatsBase

spglib = pyimport("spglib")

# lattice = [[4428.000519758899, 2484.072918780167, 0.0],
#            [-4365.2705123956675, 2592.7244786925335, 0.0],
#            [0.0, 0.0, 10000.0]]

lattice = [ [sqrt(3)/2, 1/2, 0.0],
            [-sqrt(3)/2, 1/2, 0.0],
            [0.0, 0.0, 10000]]

a1_norm = LinearAlgebra.norm(lattice[1, :])
a2_norm = LinearAlgebra.norm(lattice[2, :])
a3_norm = LinearAlgebra.norm(lattice[3, :])

angle = pi/3
d = sqrt((a1_norm^2)/(2.0*(1.0-cos(2.0*angle))))
d1 = [d*cos(angle/2.0), d*sin(angle/2.0), 0.0]

positions = [[0.0, 0.0, 0.0],
             d1]

numbers = [1, 2]
spins = [1.0, 1.0]

graphene = (lattice, positions, numbers, spins)

sym_data = spglib.get_magnetic_symmetry_dataset(graphene)

println("Group number:           ", sym_data.uni_number)
println("Hall number:            ", sym_data.hall_number)
println("Symmetry operations:    ", sym_data.n_operations)

ax1 = subplot(111, projection="3d")
ax1.scatter(lattice[1][1], lattice[1][2], lattice[1][3])
ax1.scatter(lattice[2][1], lattice[2][2], lattice[2][3])
ax1.scatter(lattice[3][1], lattice[3][2], lattice[3][3])
ax1.scatter(d1[1], d1[2], d1[3])
ax1.scatter(0.0, 0.0, 0.0)

show()

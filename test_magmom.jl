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

# PLOTTING ALL TYPES OF POINTS
latAA = read_lattice_3d("data/0.0191435_0.001/latticeAA.dat")
latBA = read_lattice_3d("data/0.0191435_0.001/latticeBA.dat")
latAB = read_lattice_3d("data/0.0191435_0.001/latticeAB.dat")
latBB = read_lattice_3d("data/0.0191435_0.001/latticeBB.dat")

bb = ase.Atoms()
for pos in eachrow(latBB)
    aux_atom = ase.Atom("C", (pos[1], pos[2], pos[3]))
    bb.append(aux_atom)
end

# ax1 = subplot(111, aspect=1)
# ax1.scatter(latAA[:,1], latAA[:,2], s=2, color="blue")
# ax1.scatter(latBA[:,1], latBA[:,2], s=2, color="green")
# ax1.scatter(latAB[:,1], latAB[:,2], s=2, color="orange")
# ax1.scatter(latBB[:,1], latBB[:,2], s=2, color="red")
#
# ax1.set_xlim([-9000, 9000])
# ax1.set_ylim([-9000, 9000])
# ax1.set_aspect("equal")
#
# show()

include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot


radius = 10000.0
latA1 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeAA.dat", radius+2000.0, radius-2000.0))
latB1 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeBA.dat", radius+2000.0, radius-2000.0))
latA2 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeAB.dat", radius+2000.0, radius-2000.0))
latB2 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeBB.dat", radius+2000.0, radius-2000.0))

ax1 = subplot(111)
ax1.scatter(latA1[1,:], latA1[2,:], s=50, color="green")
ax1.scatter(latB1[1,:], latB1[2,:], s=50, color="blue")
ax1.scatter(latA2[1,:], latA2[2,:], s=50, color="red")
ax1.scatter(latB2[1,:], latB2[2,:], s=50, color="orange")
ax1.set_aspect("equal")

write_lattice(transpose(latA1), "data/0.0191435_bernal/latticeA1_slim.dat")
write_lattice(transpose(latB1), "data/0.0191435_bernal/latticeB1_slim.dat")
write_lattice(transpose(latA2), "data/0.0191435_bernal/latticeA2_slim.dat")
write_lattice(transpose(latB2), "data/0.0191435_bernal/latticeB2_slim.dat")

show()

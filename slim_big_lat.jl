include("hex_utils.jl")
using .HexUtils


radius = 10000.0
latA1 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeA1_out.dat", radius+300.0, radius-300.0))
latB1 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeB1_out.dat", radius+300.0, radius-300.0))
# latA2 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeA2_out.dat", radius+300.0, radius-300.0))
# latB2 = transpose(read_lattice_3d("data/0.0191435_bernal/latticeB2_out.dat", radius+300.0, radius-300.0))

write_lattice(transpose(latA1), "data/0.0191435_bernal/latticeA1_slim.dat")
write_lattice(transpose(latB1), "data/0.0191435_bernal/latticeB1_slim.dat")
# write_lattice(transpose(latA2), "data/0.0191435_bernal/latticeA2_slim.dat")
# write_lattice(transpose(latB2), "data/0.0191435_bernal/latticeB2_slim.dat")

show()

include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot

# PLOTTING ALL TYPES OF POINTS
latA1 = read_lattice_3d("data/0.0191435_bernal/latticeA1_slim.dat")# , norm_a1+200.0)
latB1 = read_lattice_3d("data/0.0191435_bernal/latticeB1_slim.dat")# , norm_a1+200.0)
# latA2 = read_lattice_3d("data/0.0191435_bernal/latticeA2_slim.dat")# , norm_a1+200.0)
# latB2 = read_lattice_3d("data/0.0191435_bernal/latticeB2_slim.dat")# , norm_a1+200.0)

ax1 = subplot(111, aspect=1)
ax1.scatter(latA1[:,1], latA1[:,2], s=10, color="blue")
ax1.scatter(latB1[:,1], latB1[:,2], s=10, color="green")
# ax1.scatter(latA2[:,1], latA2[:,2], s=10, color="orange")
# ax1.scatter(latB2[:,1], latB2[:,2], s=10, color="red")

# for (i, pt) in enumerate(eachrow(latAB))
#     ax1.annotate(chop(string(latAB[i,1]), tail=8)*", "*chop(string(latAB[i,2]), tail=8), (latAB[i,1], latAB[i,2]))
# end

ax1.set_xlim([-10000, 10000])
ax1.set_ylim([-10000, 10000])
ax1.set_aspect("equal")

# legend(["A1", "B1", "A2", "B2"])

show()

include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using PyCall: LinearAlgebra


lattice = [[ 8742.840013282768, 5159.883256290869,     0.0],
           [-8840.009987150486, 4991.579924581863,     0.0],
           [               0.0,               0.0, 10000.0]]

radius = LinearAlgebra.norm(lattice[1, :])
println(radius)

# PLOTTING ALL TYPES OF POINTS
latA1 = read_lattice_3d("data/0.0191435_bernal/latticeA1_slim2.dat")# , norm_a1+200.0)
latB1 = read_lattice_3d("data/0.0191435_bernal/latticeB1_slim2.dat")# , norm_a1+200.0)
# latA2 = read_lattice_3d("data/0.0191435_bernal/latticeA2_slim.dat")# , norm_a1+200.0)
# latB2 = read_lattice_3d("data/0.0191435_bernal/latticeB2_slim.dat")# , norm_a1+200.0)

latAA = read_lattice_3d("data/0.0191435_bernal/latticeAA.dat", radius+20.0, radius-20.0)
latBA = read_lattice_3d("data/0.0191435_bernal/latticeBA.dat", radius+20.0, radius-20.0)
latAB = read_lattice_3d("data/0.0191435_bernal/latticeAB.dat", radius+20.0, radius-20.0)
latBB = read_lattice_3d("data/0.0191435_bernal/latticeBB.dat", radius+20.0, radius-20.0)

ax1 = subplot(111, aspect=1)
ax1.scatter(latA1[:,1], latA1[:,2], s=10, color="blue")
ax1.scatter(latB1[:,1], latB1[:,2], s=10, color="green")
# ax1.scatter(latA2[:,1], latA2[:,2], s=10, color="orange")
# ax1.scatter(latB2[:,1], latB2[:,2], s=10, color="red")

try ax1.scatter(latAA[:,1], latAA[:,2], s=20, color="black")
catch e
    println("No AA points")
end
try ax1.scatter(latBA[:,1], latBA[:,2], s=20, color="black")
catch e
    println("No BA points")
end
try ax1.scatter(latAB[:,1], latAB[:,2], s=20, color="black")
catch e
    printl("No AB points")
end
try ax1.scatter(latBB[:,1], latBB[:,2], s=20, color="black")
catch e
    println("No BB points")
end

# for (i, pt) in enumerate(eachrow(latAB))
#     ax1.annotate(chop(string(latAB[i,1]), tail=8)*", "*chop(string(latAB[i,2]), tail=8), (latAB[i,1], latAB[i,2]))
# end

ax1.set_xlim([-10000, 10000])
ax1.set_ylim([-10000, 10000])
ax1.set_aspect("equal")

# legend(["A1", "B1", "A2", "B2"])

show()

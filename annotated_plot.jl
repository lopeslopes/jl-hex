include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot


# PLOTTING ALL TYPES OF POINTS
latAA = read_lattice_3d("data/0.0191435_0.001/latticeAA.dat")
latBA = read_lattice_3d("data/0.0191435_0.001/latticeBA.dat")
latAB = read_lattice_3d("data/0.0191435_0.001/latticeAB.dat")
latBB = read_lattice_3d("data/0.0191435_0.001/latticeBB.dat")

ax1 = subplot(111, aspect=1)
ax1.scatter(latAA[:,1], latAA[:,2], s=20, color="blue")
ax1.scatter(latBA[:,1], latBA[:,2], s=20, color="green")
ax1.scatter(latAB[:,1], latAB[:,2], s=20, color="orange")
ax1.scatter(latBB[:,1], latBB[:,2], s=20, color="red")

for (i, pt) in enumerate(eachrow(latBA))
    ax1.annotate(chop(string(latBA[i,1]), tail=8)*", "*chop(string(latBA[i,2]), tail=8), (latBA[i,1], latBA[i,2]))
end

ax1.set_xlim([-9000, 9000])
ax1.set_ylim([-9000, 9000])
ax1.set_aspect("equal")

show()

# POINTS SELECTED FOR EACH TYPE:
# lattice = [[ 8742.840013282768, 5159.883256290869,      0.0],
#            [-8840.009987150486, 4991.579924581863,      0.0],
#            [               0.0,               0.0, 100000.0]]

# AA:
# p1 = [ 8742.840013282768, 5159.883256290869, 0.0],
# p2 = [-97.16997386763752, 10151.46318086847, 0.0],
# p3 = [-8840.009987150486, 4991.579924581863, 0.0],
# p4 = [-8742.840013282766,-5159.883256290868, 0.0],
# p5 = [ 97.1699738677083,-10151.463180872168, 0.0],
# p6 = [8840.009987150488, -4991.579924581928, 0.0]]

# BA:
# p1 = []
# p2 = []
# p3 = []
# p4 = []
# p5 = []
# p6 = []

# AB:
# p1 = []
# p2 = []
# p3 = []
# p4 = []
# p5 = []
# p6 = []

# BB:
# p1 = []
# p2 = []
# p3 = []
# p4 = []
# p5 = []
# p6 = []

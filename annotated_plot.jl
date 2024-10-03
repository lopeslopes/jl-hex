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
ax1.scatter(latAA[:,1], latAA[:,2], s=2, color="blue")
ax1.scatter(latBA[:,1], latBA[:,2], s=2, color="green")
ax1.scatter(latAB[:,1], latAB[:,2], s=2, color="orange")
ax1.scatter(latBB[:,1], latBB[:,2], s=2, color="red")

# for (i, pt) in enumerate(eachrow(latBA))
#     ax1.annotate(chop(string(latBA[i,1]), tail=8)*", "*chop(string(latBA[i,2]), tail=8), (latBA[i,1], latBA[i,2]))
# end

ax1.set_xlim([-9000, 9000])
ax1.set_ylim([-9000, 9000])
ax1.set_aspect("equal")

show()

# POINTS SELECTED FOR EACH TYPE:
#
# AA:
# lattice = [[ 8742.840013282768, 5159.883256290869, 0.0],
#            [-8840.009987150486, 4991.579924581863, 0.0],
#            [ 0.0,               0.0,               0.0]]
# positions = [[ 0.0,              0.0,               0.0],
#              [ 4314.83949352386, 2675.810337510702, 0.0],
#              [ 4428.00051975889, 2484.072918780167, 0.0],
#              [-4365.27051239566, 2592.724478692533, 0.0],
#              [-4474.73947475482, 2398.855445889395, 0.0]]
#
# BA:
# lattice = []
# position = []
#
# AB:
# lattice = []
# positions = []
#
# BB:
# lattice = []
# positions = []

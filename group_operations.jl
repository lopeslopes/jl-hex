include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Base.Threads


tol = 5.0e-3
name = @sprintf("%6.4f", tol)
origin = [0.0, 0.0]

latAA = read_lattice("data/0.0191434/"*name*"_AA.dat")
latBA = read_lattice("data/0.0191434/"*name*"_BA.dat")
latAB = read_lattice("data/0.0191434/"*name*"_AB.dat")
latBB = read_lattice("data/0.0191434/"*name*"_BB.dat")

treeAB = KDTree(transpose(latAB))
ind, dist = knn(treeAB, origin, 10)
point = latAB[ind[10],:]
println(point)

slope = point[2]/point[1]
x_vals = -12500:12500
y_vals = slope*x_vals

ax1 = subplot(111, aspect=1)
ax1.scatter(latAA[:,1], latAA[:,2], color="green", s=8)
ax1.scatter(latBA[:,1], latBA[:,2], color="red", s=8)
ax1.scatter(latAB[:,1], latAB[:,2], color="orange", s=8)
ax1.scatter(latBB[:,1], latBB[:,2], color="blue", s=8)

#ax1.quiver(0.0, 0.0, point[1], point[2], angles="xy", scale_units="xy", scale=1)
ax1.plot(x_vals, y_vals)

ax1.set_xlim([-12250, 12250])
ax1.set_ylim([-12250, 12250])
legend(["AA", "BA", "AB", "BB"])

show()

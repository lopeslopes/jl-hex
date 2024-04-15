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

#latAA = read_lattice("data/0.0191434/"*name*"_AA.dat")
latBA = read_lattice("data/0.0191434/"*name*"_BA.dat")
latAB = read_lattice("data/0.0191434/"*name*"_AB.dat")
#latBB = read_lattice("data/0.0191434/"*name*"_BB.dat")

new_points = []
for i in 1:size(latAB,1)
    point = latAB[i,:]
    np = sym_op_d6h(point, 10)
    push!(new_points, np)
end
new_lat = transpose(hcat(new_points...))

ax1 = subplot(111, aspect=1)
#ax1.scatter(latAA[:,1], latAA[:,2], color="green", s=8)
ax1.scatter(latBA[:,1], latBA[:,2], color="blue", s=8)
ax1.scatter(latAB[:,1], latAB[:,2], color="green", s=8)
#ax1.scatter(latBB[:,1], latBB[:,2], color="green", s=8)

ax1.scatter(new_lat[:,1], new_lat[:,2], color="red", s=1)

show()

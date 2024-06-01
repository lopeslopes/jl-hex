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
ind, dist = knn(treeAB, origin, 1000)
lenBB = size(latBB)[1]

point_AB = [0.0, 0.0]
point_BB = [0.0, 0.0]
slope_final = 0.0
for i in 1:1000
    point = latAB[ind[i],:]
    slope = point[2]/point[1]
    for j in 1:lenBB
        y_slope = latBB[j,1]*slope
        diff = abs(latBB[j,2]-y_slope)
        if diff<0.1
            println("AB point: ")
            println("    index:  ", i)
            println("    coords: ", point)
            println("BB point: ")
            println("    index:  ", j)
            println("    coords: ", latBB[j,:])
            println("    diff:   ", diff)
            global point_AB = point
            global point_BB = latBB[j,:]
            global slope_final = point_AB[2]/point_AB[1]
            println("Slope of line: ", slope_final)
        end
    end
end

x_vals = -12500:12500
y_vals = slope_final*x_vals

ax1 = subplot(111, aspect=1)
ax1.scatter(latAA[:,1], latAA[:,2], color=("green", 0.5), s=8)
ax1.scatter(latBA[:,1], latBA[:,2], color=("red", 0.5), s=8)
ax1.scatter(latAB[:,1], latAB[:,2], color=("orange", 0.5), s=8)
ax1.scatter(latBB[:,1], latBB[:,2], color=("blue", 0.5), s=8)

# SELECTED AB and BB points
ax1.scatter(0.0, 0.0, color="green", s=50)
ax1.scatter(point_AB[1], point_AB[2], color="magenta", s=50)
ax1.scatter(point_BB[1], point_BB[2], color="purple", s=50)

ax1.plot(x_vals, y_vals)

ax1.set_xlim([-12250, 12250])
ax1.set_ylim([-12250, 12250])
legend(["AA", "BA", "AB", "BB", "selected AA", "selected AB", "selected BB"])

show()

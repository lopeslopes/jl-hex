include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using PyCall: LinearAlgebra
using NearestNeighbors
using DelaunayTriangulation


latAA = read_lattice("data/0.3802512/latticeAA.dat")
latBA = read_lattice("data/0.3802512/latticeBA.dat")
latAB = read_lattice("data/0.3802512/latticeAB.dat")
latBB = read_lattice("data/0.3802512/latticeBB.dat")


tree = KDTree(transpose(latAB))
ind_nbors, dist_nbors = knn(tree, [0.0, 0.0], 7)
neighbors = [latAB[i,:] for i in ind_nbors]

tri = triangulate(neighbors)
tess = voronoi(tri)

poly = [[pt[1], pt[2]] for pt in tess.polygon_points]
tree = KDTree(hcat(poly...))
poly_ord = []
used_ind = []
flag = false
i = 1
while (flag == false)
    push!(poly_ord, poly[i][:])
    push!(used_ind, i)
    ind, dist = knn(tree, [poly[i][1], poly[i][2]], 3)
    deleteat!(ind, argmin(dist))
    deleteat!(dist, argmin(dist))
    if (size(used_ind) == size(poly))
        global flag = true
    end
    if (ind[1] in used_ind)
        global i = ind[2]
    else
        global i = ind[1]
    end
end

ax1 = subplot(111, aspect=1)

poly_ord = hcat(poly_ord...)

x = poly_ord[1,:]
y = poly_ord[2,:]
println(x)
println(y)

poly_ord = transpose(poly_ord)
ax1.plot(poly_ord)



# try ax1.scatter(latAA[:,1], latAA[:,2], s=20, color="blue")
# catch e
#     println("No AA points")
# end
# try ax1.scatter(latBA[:,1], latBA[:,2], s=20, color="orange")
# catch e
#     println("No BA points")
# end
# try ax1.scatter(latAB[:,1], latAB[:,2], s=20, color="purple")
# catch e
#     printl("No AB points")
# end
# try ax1.scatter(latBB[:,1], latBB[:,2], s=20, color="magenta")
# catch e
#     println("No BB points")
# end

ax1.scatter(x, y)

ax1.set_xlim([-10, 10])
ax1.set_ylim([-10, 10])
ax1.set_aspect("equal")

show()

include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using PyCall: LinearAlgebra
using NearestNeighbors
using DelaunayTriangulation


latAA = read_lattice("data/0.0189290/latticeAA.dat")
latBA = read_lattice("data/0.0189290/latticeBA.dat")
latAB = read_lattice("data/0.0189290/latticeAB.dat")
latBB = read_lattice("data/0.0189290/latticeBB.dat")


origin = [0.0, 0.0]
tree = KDTree(transpose(latAB))
ind_nbors, dist_nbors = knn(tree, origin, 7)
neighbors = [latAB[i,:] for i in ind_nbors]

ax1 = subplot(111, aspect=1)
legend = []

## Wigner Seitz cell
# tri = triangulate(neighbors)
# tess = voronoi(tri)
#
# poly = [[pt[1], pt[2]] for pt in tess.polygon_points]
# tree = KDTree(hcat(poly...))
# poly_ord = []
# used_ind = []
# flag = false
# i = 1
# while (flag == false)
#     push!(poly_ord, poly[i][:])
#     push!(used_ind, i)
#     ind, dist = knn(tree, [poly[i][1], poly[i][2]], 3)
#     deleteat!(ind, argmin(dist))
#     deleteat!(dist, argmin(dist))
#     if (size(used_ind) == size(poly))
#         global flag = true
#     end
#     if (ind[1] in used_ind)
#         global i = ind[2]
#     else
#         global i = ind[1]
#     end
# end
# push!(poly_ord, poly_ord[1][:])
# poly_ord = hcat(poly_ord...)
#
# x = poly_ord[1,:]
# y = poly_ord[2,:]
#
# ax1.plot(x, y)
# push!(legend, "Wigner-Seitz cell")

# area_ws_cell = DelaunayTriangulation.get_largest_area(DelaunayTriangulation.statistics(tri))
# n_poly = DelaunayTriangulation.num_polygons(tess)
# area_ws_cell = 0.0
# for i in 1:n_poly
#     area = DelaunayTriangulation.polygon_features(tess, i)[1]
#     if (area != Inf)
#         global area_ws_cell = area
#     end
# end
# println("Area (wigner seitz): ", area_ws_cell)


## lattice vectors cell
a1 = neighbors[1][:]
a2 = copy(a1)
rotate_point!(a2, (2.0/3.0)*pi, origin)
a3 = a2 + a1
poly_ord = [origin,
            a1,
            a3,
            a2,
            origin]
poly_ord = hcat(poly_ord...)
x = poly_ord[1,:]
y = poly_ord[2,:]

ax1.plot(x, y)
push!(legend, "Lattice vectors cell")

a1_3d = [a1[1], a1[2], 0.0]
a2_3d = [a2[1], a2[2], 0.0]
area_lat_cell = LinearAlgebra.norm(LinearAlgebra.cross(a1_3d, a2_3d))
println("Area (lattice vecs): ", area_lat_cell)

println("Lat. vector obtained: ", LinearAlgebra.norm(a1))

# superlattice vectors determination by article formula
ang_aux = (60.0/180.0) * π
a = 2.46
v1 = [a, 0.0]
v2 = [a*cos(ang_aux), a*sin(ang_aux)]

p = 1
q = 61
angle = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
println("Angle in radians: ", angle)
# println("Angle in degrees: ", (angle * 180) / pi)

delta = 3/gcd(p, 3)
gamma = gcd(p+3*q, p-3*q)
println("Delta: ",delta)
println("Gamma: ", gamma)

if (delta == 1.0)
    t1 = (1/gamma)*(p+3*q)*v1 + (1/gamma)*(p-3*q)*v2
    t2 = -(1/gamma)*(2*p)*v1 + (1/gamma)*(p+3*q)*v2
elseif(delta == 3.0)
    t1 = (1/gamma)*(p+q)*v1 - (1/gamma)*(2*q)*v2
    t2 = -(1/gamma)*(p-q)*v1 + (1/gamma)*(p+q)*v2 
end

println("|t1| = ", LinearAlgebra.norm(t1))
println("|t2| = ", LinearAlgebra.norm(t2))

moire_period = a/(2*sin(angle/2))
println("D = ", moire_period)

# try 
#     ax1.scatter(latAA[:,1], latAA[:,2], s=10, color="blue")
#     push!(legend, "AA points")
# catch e
#     println("No AA points")
# end
# try 
#     ax1.scatter(latBA[:,1], latBA[:,2], s=10, color="orange")
#     push!(legend, "BA points")
# catch e
#     println("No BA points")
# end
# try 
#     ax1.scatter(latAB[:,1], latAB[:,2], s=70, color="purple")
#     push!(legend, "AB points")
# catch e
#     println("No AB points")
# end
# try 
#     ax1.scatter(latBB[:,1], latBB[:,2], s=10, color="magenta")
#     push!(legend, "BB points")
# catch e
#     println("No BB points")
# end
#
# ax1.set_xlim([-400, 400])
# ax1.set_ylim([-400, 400])
# ax1.set_aspect("equal")

# ax1.legend(legend)

# show()

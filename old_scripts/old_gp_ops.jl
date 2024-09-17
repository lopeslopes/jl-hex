include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Spglib

tol = 5.0e-3
name = @sprintf("%6.4f", tol)
origin = [0.0, 0.0]

latAA = read_lattice("data/0.0191434/"*name*"_AA.dat")
latBA = read_lattice("data/0.0191434/"*name*"_BA.dat")
latAB = read_lattice("data/0.0191434/"*name*"_AB.dat")
latBB = read_lattice("data/0.0191434/"*name*"_BB.dat")

treeAA = KDTree(transpose(latAA))

ind, dist = knn(treeAA, origin, 10)
point = latAA[ind[1],:]

sym_points = []
for i in 1:24
    p_aux = sym_op_d6h(point, i)
    push!(sym_points, p_aux)
end
sym_points = transpose(hcat(sym_points...))

treeBA = KDTree(transpose(latBA))
treeAB = KDTree(transpose(latAB))
treeBB = KDTree(transpose(latBB))

cAA = 0
cBA = 0
cAB = 0
cBB = 0
tol_sym = 0.00001
for j in 1:24
    sp = sym_points[j,:]
    ind1, dist1 = knn(treeAA, sp, 1)
    if dist1[1]>tol_sym
        global cAA += 1
    end
    ind1, dist1 = knn(treeBA, sp, 1)
    if dist1[1]>tol_sym
        global cBA += 1
    end
    ind1, dist1 = knn(treeAB, sp, 1)
    if dist1[1]>tol_sym
        global cAB += 1
    end
    ind1, dist1 = knn(treeBB, sp, 1)
    if dist1[1]>tol_sym
        global cBB += 1
    end
end
println("cAA = ", cAA)
println("cBA = ", cBA)
println("cAB = ", cAB)
println("cBB = ", cBB)

#ax1 = subplot(111, aspect="equal")
#ax1.scatter(sym_points[:,1], sym_points[:,2], s=20, color="c")
#
#show()

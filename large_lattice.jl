include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Base.Threads


# INITIAL DEFINITIONS
n = 200000000
a = 2.46
origin1 = [0.0, 0.0]
origin2 = [0.0, 0.0]

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)
latA2 = zeros(n ÷ 2, 2)
latB2 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)
HexUtils.create_honeycomb_lattice!(latA2, latB2, a, false)

angle = 1.91517084211538477391489120822979632e-2
tol = 5.0e-3

println("Tolerance:        ", tol)
println("Angle in radians: ", angle)
println("Angle in degrees: ", (angle * 180) / pi)

# ROTATE SECOND LATTICE BY THE ANGLE
rotate_lattice!(latA2, angle, origin2)
rotate_lattice!(latB2, angle, origin2)

# TEST SECTION: TREES
treeA1 = KDTree(transpose(latA1))
treeB1 = KDTree(transpose(latB1))

AA = []
BA = []
AB = []
BB = []

for i in 1:div(n,2)
    indAA, distAA = knn(treeA1, latA2[i,:], 1)
    indBA, distBA = knn(treeB1, latA2[i,:], 1)
    indAB, distAB = knn(treeA1, latB2[i,:], 1)
    indBB, distBB = knn(treeB1, latB2[i,:], 1)
    if distAA[1] < tol
        push!(AA, latA2[i,:])
    end
    if distBA[1] < tol
        push!(BA, latA2[i,:])
    end
    if distAB[1] < tol
        push!(AB, latB2[i,:])
    end
    if distBB[1] < tol
        push!(BB, latB2[i,:])
    end
end

latAA = transpose(hcat(AA...))
latBA = transpose(hcat(BA...))
latAB = transpose(hcat(AB...))
latBB = transpose(hcat(BB...))

tol_str = @sprintf("%6.4f", tol)

ax1 = subplot(111,aspect=1)
ax1.scatter(latAA[:,1], latAA[:,2], s=1, color="blue")
ax1.scatter(latBA[:,1], latBA[:,2], s=1, color="green")
ax1.scatter(latAB[:,1], latAB[:,2], s=1, color="orange")
ax1.scatter(latBB[:,1], latBB[:,2], s=1, color="red")

legend(["AA", "BA", "AB", "BB"])
savefig("results/1.9151/"*tol_str*"_200M.png", format="png", dpi=550)

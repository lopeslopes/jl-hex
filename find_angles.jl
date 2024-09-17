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
origin = [0.0, 0.0]

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)
latA2 = zeros(n ÷ 2, 2)
latB2 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)
HexUtils.create_honeycomb_lattice!(latA2, latB2, a, false)

distsA1 = zeros(n ÷ 2)
distsB1 = zeros(n ÷ 2)
distsA2 = zeros(n ÷ 2)
distsB2 = zeros(n ÷ 2)
for i in 1:div(n,2)
    distsA1[i] = sqrt(latA1[i,1]^2 + latA1[i,2]^2)
    distsA2[i] = sqrt(latA2[i,1]^2 + latA2[i,2]^2)
    distsB1[i] = sqrt(latB1[i,1]^2 + latB1[i,2]^2)
    distsB2[i] = sqrt(latB2[i,1]^2 + latB2[i,2]^2)
end

# TEST SECTION: TREES
treeA1 = KDTree(transpose(distsA1), reorder=false)
treeB1 = KDTree(transpose(distsB1), reorder=false)
 
tol = 1.0e-4
BA = []
AB = []
println("Tolerance:        ", tol)
for i in 1:div(n,2)
    indBA, distBA = nn(treeB1, [distsA2[i]])
    indAB, distAB = nn(treeA1, [distsB2[i]])
    if distBA < tol
        angA = atan(latA2[i,2]/latA2[i,1])
        angB = atan(latB1[indBA,2]/latB1[indBA,1])
        angle = max(angA, angB) - min(angA, angB)
        if 0.0191340 < angle < 0.0191430
            push!(BA, angle)
        end
    end
    if distAB < tol
        angA = atan(latA1[indAB,2]/latA1[indAB,1])
        angB = atan(latB2[i,2]/latB2[i,1])
        angle = max(angA, angB) - min(angA, angB)
        if 0.0191340 < angle < 0.0191430
            push!(AB, angle)
        end
    end
end

println(sort(union(BA..., AB...)))

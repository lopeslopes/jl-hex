include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Base.Threads


# INITIAL DEFINITIONS
n = 120000000
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

# NEW ANGLES OBTAINED FOR TESTING:
# 0.019120090428629055
# 0.019124340243658833
# 0.01912540230893478
# 0.01914345108312343
# 0.019151034345471163
# 0.01915472274745922
# 0.019155348627157176
# 0.019167605306979363
#
#
# 0.019132574878803932
# 0.01913264068583198
# 0.019133982309114894
# 0.019134082825285392
# 0.019134209451742046
# 0.01914204154491761
# 0.01914213838482104
# 0.019143422530969456
# 0.01914348511903463
#
# 0.01913408282528528
# 0.01913420945174149
# 0.019134299302355773
# 0.019134412838422543
# 0.01913449363522135
# 0.019141644839413052
# 0.019141722979244657
# 0.01914183268481051
# 0.019141919423809495
# 0.019142041544964017
# 0.019142138384832696

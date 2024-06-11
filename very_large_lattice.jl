include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf


# INITIAL DEFINITIONS
n = 360000000
a = 2.46
hex_center_pivot = false
AB_stacking = false

# STACKING AND ORIGIN DEFINITION
if AB_stacking
    println("Stacking mode: AB (Bernal stacking)")
else
    println("Stacking mode: AA (no displacement)")
end

if hex_center_pivot
    println("Pivot point: empty center of hexagonal cell")
    d = sqrt((a^2) / (2 * (1 - cos(2 * pi / 3))))
    d1 = [d * cos(pi / 6), d * sin(pi / 6)]
    origin1 = d1 - [0.0, 0.0]
    origin2 = d1
else
    println("Pivot point: node at origin")
    origin1 = [0.0, 0.0]
    origin2 = [0.0, 0.0]
end

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)
latA2 = zeros(n ÷ 2, 2)
latB2 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)
HexUtils.create_honeycomb_lattice!(latA2, latB2, a, AB_stacking)

angle = 0.01914345108312343
println("Angle in radians: ", angle)
println("Angle in degrees: ", (angle * 180) / pi)

# ROTATE SECOND LATTICE BY THE ANGLE
rotate_lattice!(latA2, angle, origin2)
rotate_lattice!(latB2, angle, origin2)

tol = 5.0e-3
println("Tolerance:        ", tol)
name = @sprintf("%6.4f", tol)

# TEST SECTION: TREES

AA = []
AB = []
tree = KDTree(transpose(latA1))
for i in 1:div(n,2)
    indAA, distAA = knn(tree, latA2[i,:], 1)
    indAB, distAB = knn(tree, latB2[i,:], 1)
    if distAA[1] < tol
        push!(AA, latA2[i,:])
    end
    if distAB[1] < tol
        push!(AB, latB2[i,:])
    end
end
latAA = transpose(hcat(AA...))
latAB = transpose(hcat(AB...))
write_lattice(latAA, "data/0.0191434/300M_"*name*"_AA.dat")
write_lattice(latAB, "data/0.0191434/300M_"*name*"_AB.dat")

BA = []
BB = []
tree = KDTree(transpose(latB1))
for i in 1:div(n,2)
    indBA, distBA = knn(tree, latA2[i,:], 1)
    indBB, distBB = knn(tree, latB2[i,:], 1)
    if distBA[1] < tol
        push!(BA, latA2[i,:])
    end
    if distBB[1] < tol
        push!(BB, latB2[i,:])
    end
end
latBA = transpose(hcat(BA...))
latBB = transpose(hcat(BB...))
write_lattice(latBA, "data/0.0191434/300M_"*name*"_BA.dat")
write_lattice(latBB, "data/0.0191434/300M_"*name*"_BB.dat")

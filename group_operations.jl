include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf

ase = pyimport("ase")
spglib = pyimport("spglib")

tol = 5.0e-3
name = @sprintf("%6.4f", tol)
origin = [0.0, 0.0]

# INITIAL DEFINITIONS
n = 500000
a = 2.46
hex_center_pivot = false
AB_stacking = true

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

tol = 5.0e-3
name = @sprintf("%6.4f", tol)

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating bottom lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)

println("Creating top lattices...")
latA2 = zeros(n ÷ 2, 2)
latB2 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA2[:,1:2], latB2[:,1:2], a, AB_stacking)

angle = 0.01914345108312343
println("Angle in radians: ", angle)
println("Angle in degrees: ", (angle * 180) / pi)

# ROTATE SECOND LATTICE BY THE ANGLE
rotate_lattice!(latA2[:,1:2], angle, origin2)
rotate_lattice!(latB2[:,1:2], angle, origin2)

tol = 1.0e-3
println("Tolerance:        ", tol)
name = @sprintf("%6.4f", tol)

AA = []
BA = []
AB = []
BB = []

treeA1 = KDTree(transpose(latA1))
treeB1 = KDTree(transpose(latB1))

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

ang_name = @sprintf("%9.7f", angle)
mkdir("data/"*ang_name)

try
    write_lattice(latAA, "data/"*ang_name*"/"*name*"_AA.dat")
catch e
    println("AA lattice is empty!")
end
try
    write_lattice(latBA, "data/"*ang_name*"/"*name*"_BA.dat")
catch e
    println("BA lattice is empty!")
end
try
    write_lattice(latAB, "data/"*ang_name*"/"*name*"_AB.dat")
catch e
    println("AB lattice is empty!")
end
try
    write_lattice(latBB, "data/"*ang_name*"/"*name*"_BB.dat")
catch e
    println("BB lattice is empty!")
end

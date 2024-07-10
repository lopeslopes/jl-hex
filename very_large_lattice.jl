include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf


# INITIAL DEFINITIONS
n = 200000000
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

tol = 5.0e-3
name = @sprintf("%6.4f", tol)

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating bottom lattices...")
latA = zeros(n ÷ 2, 2)
latB = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA, latB, a, false)

write_lattice(latA, "data/0.0191434/1B_"*name*"_A1.dat")
write_lattice(latB, "data/0.0191434/1B_"*name*"_B1.dat")

println("Creating top lattices...")
angle = 0.01914345108312343
println("Angle in radians: ", angle)
println("Angle in degrees: ", (angle * 180) / pi)

# ROTATE SECOND LATTICE BY THE ANGLE
rotate_lattice!(latA, angle, origin2)
rotate_lattice!(latB, angle, origin2)

write_lattice(latA, "data/0.0191434/1B_"*name*"_A2.dat")
write_lattice(latB, "data/0.0191434/1B_"*name*"_B2.dat")

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
n = 200
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

HexUtils.create_honeycomb_lattice_fractional!(latA1, latB1, a, false)

#println("Creating top lattices...")
#latA2 = zeros(n ÷ 2, 2)
#latB2 = zeros(n ÷ 2, 2)
#
#HexUtils.create_honeycomb_lattice!(latA2[:,1:2], latB2[:,1:2], a, AB_stacking)
#
#angle = 0.01914345108312343
#println("Angle in radians: ", angle)
#println("Angle in degrees: ", (angle * 180) / pi)
#
## ROTATE SECOND LATTICE BY THE ANGLE
#rotate_lattice!(latA2[:,1:2], angle, origin2)
#rotate_lattice!(latB2[:,1:2], angle, origin2)

lenA1 = size(latA1)[1]
typeA1 = ones(Int64, lenA1)
println("Lattices created, starting symmetry calculation")

lat_angle = pi/3.0
a1 = [a, 0.0, 0.0]
a2 = [a*cos(lat_angle), a*sin(lat_angle), 0.0]
a3 = [0.0, 0.0, 1.0]

# CONVERTING CARTESIAN TO FRACTIONAL COORDINATES
#latA1_frac = cartesian_to_fractional(hcat(a1,a2,a3)..., latA1)
latA1_3d = zeros(n ÷ 2, 3)
for i=1:lenA1
    latA1_3d[i,:] = [latA1[i,1], latA1[i,2], 0.0]
end

test_lat = spglib.Lattice([a1,a2,a3])
test_cell = spglib.Cell(test_lat, latA1_3d, typeA1)
println(test_lat)
println(test_cell)
test = spglib.get_symmetry(test_cell)

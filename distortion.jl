include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Base.Threads


# INITIAL DEFINITIONS
n = 10000000
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

angle = 1.90495757117080288768798873954979561e-2
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
treeA2 = KDTree(transpose(latA2))
treeB2 = KDTree(transpose(latB2))

# DISTORTION OF ALL LAT2 BY THE SAME DIFF VECTOR
# test_point = [ 33.20919997773677, 452.36422740621947]
# test_point = [-41.82071381645537, 451.64503616874754]
# test_point = [ -2.46429740059443, 258.49122158850423]

# ind_test, dist_test = knn(treeB1, test_point, 1)
# foundB1 = latB1[ind_test[1],:]
# ind_test, dist_test = knn(treeA2, test_point, 1)
# foundA2 = latA2[ind_test[1],:]
# 
# diff = foundB1 - foundA2
diff = [0.0, 0.0]

latA2[:,1] = latA2[:,1] .+ diff[1]
latA2[:,2] = latA2[:,2] .+ diff[2]
latB2[:,1] = latB2[:,1] .+ diff[1]
latB2[:,2] = latB2[:,2] .+ diff[2]

# NEW OVERLAP TEST AFTER DISTORTION
tol = 5.0e-3
println("New overlap test after distortion starting...")
println("Tolerance:        ", tol)

treeA2 = KDTree(transpose(latA2))
treeB2 = KDTree(transpose(latB2))

AA = []
BA = []
AB = []
BB = []

#Threads.nthreads() = 6
#@threads for i in 1:div(n,2)
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

latAA = hcat(AA...)
latBA = hcat(BA...)
latAB = hcat(AB...)
latBB = hcat(BB...)

# PLOT
ax1 = subplot(111,aspect=1)
ax1.scatter(latAA[1,:], latAA[2,:], color="blue")
ax1.scatter(latBA[1,:], latBA[2,:], color="green")
ax1.scatter(latAB[1,:], latAB[2,:], color="orange")
ax1.scatter(latBB[1,:], latBB[2,:], color="red")

# ax1.quiver(0.0, 0.0, foundA2[1], foundA2[2], angles="xy", scale_units="xy", scale=1)

ax1.set_xlim([-2400, 2400])
ax1.set_ylim([-2400, 2400])
legend(["AA", "BA", "AB", "BB"])
show()

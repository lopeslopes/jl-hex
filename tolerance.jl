include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf


# INITIAL DEFINITIONS
n = 150000000
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

# TEST ANGLES OBTAINED USING calc_angle.f90
# angle = 1.90264929971534712693433092186275811e-2
# angle = 1.90403102830304951869456657448196344e-2
# angle = 1.90495757117080288768798873954979561e-2
# angle = 1.91235783233762423239743767417441124e-2
# angle = 1.91517084211538477391489120822979632e-2
# angle = 1.91756277374731560763308567863988621e-2

angle = 1.91235783233762423239743767417441124e-2
println("Angle in radians: ", angle)
println("Angle in degrees: ", (angle * 180) / pi)

# ROTATE SECOND LATTICE BY THE ANGLE
rotate_lattice!(latA2, angle, origin2)
rotate_lattice!(latB2, angle, origin2)

# write_lattice(latA1, "data/1.9175/latticeA1.dat")
# write_lattice(latB1, "data/1.9175/latticeB1.dat")
# write_lattice(latA2, "data/1.9175/latticeA2.dat")
# write_lattice(latB2, "data/1.9175/latticeB2.dat")

# TEST SECTION: TREES
treeA1 = KDTree(transpose(latA1))
treeB1 = KDTree(transpose(latB1))
# treeA2 = KDTree(transpose(latA2))
# treeB2 = KDTree(transpose(latB2))

# TOLERANCE VARIATION
tols = [5.0e-3, 4.0e-3, 3.0e-3, 2.0e-3, 1.0e-3]
for tol in tols
    println("Tolerance:        ", tol)

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

    latAA = transpose(hcat(AA...))
    latBA = transpose(hcat(BA...))
    latAB = transpose(hcat(AB...))
    latBB = transpose(hcat(BB...))

    name = @sprintf("%6.4f", tol)
    # write_lattice(latAA, "data/1.9175/latAA_"*name*".dat")
    # write_lattice(latBA, "data/1.9175/latBA_"*name*".dat")
    # write_lattice(latAB, "data/1.9175/latAB_"*name*".dat")
    # write_lattice(latBB, "data/1.9175/latBB_"*name*".dat")

    # PLOT
    ax1 = subplot(111,aspect=1)
    ax1.scatter(latAA[:,1], latAA[:,2], s=1, color="blue")
    ax1.scatter(latBA[:,1], latBA[:,2], s=1, color="green")
    ax1.scatter(latAB[:,1], latAB[:,2], s=1, color="orange")
    ax1.scatter(latBB[:,1], latBB[:,2], s=1, color="red")

    ax1.set_xlim([-12250, 12250])
    ax1.set_ylim([-12250, 12250])
    legend(["AA", "BA", "AB", "BB"])
    savefig("results/1.9123/150M/"*name*"_150M.png", format="png", dpi=550)
    clf()
end

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
origin = [0.0, 0.0]

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)

tol = 4.0e-3
println("Tolerance:        ", tol)

treeA1 = KDTree(transpose(latA1))
treeB1 = KDTree(transpose(latB1))

# ANGLE VARIATION
angles = [0.01913408282528528, 0.01913420945174149, 0.019134299302355773, 0.019134412838422543, 0.01913449363522135, 0.019141644839413052, 0.019141722979244657, 0.01914183268481051, 0.019141919423809495, 0.019142041544964017, 0.019142138384832696]

for angle in angles
    println("Angle in radians: ", angle)
    println("Angle in degrees: ", (angle * 180) / pi)

    latA2 = zeros(n ÷ 2, 2)
    latB2 = zeros(n ÷ 2, 2)
    HexUtils.create_honeycomb_lattice!(latA2, latB2, a, AB_stacking)
    rotate_lattice!(latA2, angle, origin)
    rotate_lattice!(latB2, angle, origin)

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

    name = @sprintf("%19.17f", angle)

    # PLOT
    ax1 = subplot(111,aspect=1)
    ax1.scatter(latAA[:,1], latAA[:,2], s=1, color="blue")
    ax1.scatter(latBA[:,1], latBA[:,2], s=1, color="green")
    ax1.scatter(latAB[:,1], latAB[:,2], s=1, color="orange")
    ax1.scatter(latBB[:,1], latBB[:,2], s=1, color="red")

    ax1.set_xlim([-12250, 12250])
    ax1.set_ylim([-12250, 12250])
    legend(["AA", "BA", "AB", "BB"])
    savefig("results/angle_test/"*name*".png", format="png", dpi=550)
    clf()
end

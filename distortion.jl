include("hex_utils.jl")
using .HexUtils
using NearestNeighbors
using Printf


base_data_path = "data"
entries = readdir(base_data_path; join=true)
dataset_dirs = sort(filter(isdir, entries))
n_datasets = length(dataset_dirs)

path = last(dataset_dirs)
println(path)
println(" ")

latA1 = transpose(read_lattice_3d(path*"/latticeA1.dat"))
latB1 = transpose(read_lattice_3d(path*"/latticeB1.dat"))
latA2 = transpose(read_lattice_3d(path*"/latticeA2.dat"))
latB2 = transpose(read_lattice_3d(path*"/latticeB2.dat"))

treeA1 = KDTree(latA1)
treeB1 = KDTree(latB1)
treeA2 = KDTree(latA2)
treeB2 = KDTree(latB2)

latA1 = transpose(latA1)
latB1 = transpose(latB1)
latA2 = transpose(latA2)
latB2 = transpose(latB2)

AA = []
BA = []
AB = []
BB = []

n = minimum([size(latA1)[1], size(latB1)[1], size(latA2)[1], size(latB2)[1]])

smallest_AA_separation = 10000.0
AA_vector_separation = [0.0, 0.0]

smallest_BA_separation = 10000.0
BA_vector_separation = [0.0, 0.0]

smallest_AB_separation = 10000.0
AB_vector_separation = [0.0, 0.0]

smallest_BB_separation = 10000.0
BB_vector_separation = [0.0, 0.0]

for i in 1:div(n,2)
    indAA, distAA = knn(treeA1, latA2[i,:], 1)
    if distAA[1] < smallest_AA_separation
        global smallest_AA_separation = distAA[1]
        global AA_vector_separation = transpose(latA1[indAA,1:3]) - latA2[i,:]
        global AA_vector_separation = AA_vector_separation[1:3]
    end

    indBA, distBA = knn(treeB1, latA2[i,:], 1)
    if distBA[1] < smallest_BA_separation
        global smallest_BA_separation = distBA[1]
        global BA_vector_separation = transpose(latB1[indBA,1:3]) - latA2[i,:]
        global BA_vector_separation = BA_vector_separation[1:3]
    end

    indAB, distAB = knn(treeA1, latB2[i,:], 1)
    if distAB[1] < smallest_AB_separation
        global smallest_AB_separation = distAB[1]
        global AB_vector_separation = transpose(latA1[indAB,1:3]) - latB2[i,:]
        global AB_vector_separation = AB_vector_separation[1:3]
    end

    indBB, distBB = knn(treeB1, latB2[i,:], 1)
    if distBB[1] < smallest_BB_separation
        global smallest_BB_separation = distBB[1]
        global BB_vector_separation = transpose(latB1[indBB,1:3]) - latB2[i,:]
        global BB_vector_separation = BB_vector_separation[1:3]
    end
end

println("Smallest AA separation: ")
println(smallest_AA_separation)
println(AA_vector_separation)
println(" ")

println("Smallest BA separation: ")
println(smallest_BA_separation)
println(BA_vector_separation)
println(" ")

println("Smallest AB separation: ")
println(smallest_AB_separation)
println(AB_vector_separation)
println(" ")

println("Smallest BB separation: ")
println(smallest_BB_separation)
println(BB_vector_separation)
println(" ")


# latAA = transpose(hcat(AA...))
# latBA = transpose(hcat(BA...))
# latAB = transpose(hcat(AB...))
# latBB = transpose(hcat(BB...))

# if isempty(latAA)
#     println("AA: empty")
# else
#     println("AA: ", size(latAA)[1])
# end
#
# if isempty(latBA)
#     println("BA: empty")
# else
#     println("BA: ", size(latBA)[1])
# end
#
# if isempty(latAB)
#     println("AB: empty")
# else
#     println("AB: ", size(latAB)[1])
# end
#
# if isempty(latBB)
#     println("BB: empty")
# else
#     println("BB: ", size(latBB)[1])
# end
#
# # INITIAL DEFINITIONS
# n = 200000
#
# a_top = 2.46
# a1_top = [a_top, 0.0]
# a2_top = [-a_top*cos(pi/3.0), a_top*sin(pi/3.0)]
#
# AB_stacking = false
#
# rot_axis = [0.0, 0.0]
#
# # ALLOCATION OF LATTICES AND FIRST CREATION
# println("Creating lattices...")
# latA1 = zeros(n ÷ 2, 2)
# latB1 = zeros(n ÷ 2, 2)
#
# HexUtils.create_honeycomb_lattice!(latA1, latB1, a_top, a1_top, a2_top, false)
#
# max_radius = maximum(latA1) - 10.0
#
# treeA1 = KDTree(transpose(latA1))
# treeB1 = KDTree(transpose(latB1))
#
# for q in [63.0, 62.0, 61.0, 60.0, 59.0, 58.0, 57.0, 56.0]
#     p = 1.0
#     # q = 56.0
#     angle_i = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
#     angle_f = acos((3.0*((q-1)^2) - (p^2))/(3.0*((q-1)^2) + (p^2)))
#     println(angle_i)
#     println(angle_f)
#     steps = 10
#
#     ## 10 steps between each magic angle (maybe too much, but we'll see)
#     for j in 1:steps
#         latA2 = zeros(n ÷ 2, 2)
#         latB2 = zeros(n ÷ 2, 2)
#
#         a_bot = 2.46
#         a1_bot = [a_bot, 0.0]
#         a2_bot = [-a_bot*cos(pi/3.0), a_bot*sin(pi/3.0)]
#
#         HexUtils.create_honeycomb_lattice!(latA2, latB2, a_bot, a1_bot, a2_bot, AB_stacking)
#
#         angle = angle_i + j*(angle_f - angle_i)/10
#         ang_name = @sprintf("%9.7f", angle)
#         println("Angle in radians: ", angle)
#         println("Angle in degrees: ", (angle * 180) / pi)
#
#         # ROTATE SECOND LATTICE BY THE ANGLE
#         rotate_lattice!(latA2, angle, rot_axis)
#         rotate_lattice!(latB2, angle, rot_axis)
#
#         tol = 5.0e-3
#         println("Tolerance:        ", tol)
#         name = @sprintf("%6.4f", tol)
#
#         AA = []
#         BA = []
#         AB = []
#         BB = []
#
#         for i in 1:div(n,2)
#             indAA, distAA = knn(treeA1, latA2[i,:], 1)
#             indBA, distBA = knn(treeB1, latA2[i,:], 1)
#             indAB, distAB = knn(treeA1, latB2[i,:], 1)
#             indBB, distBB = knn(treeB1, latB2[i,:], 1)
#             if distAA[1] < tol
#                 push!(AA, latA2[i,:])
#             end
#             if distBA[1] < tol
#                 push!(BA, latA2[i,:])
#             end
#             if distAB[1] < tol
#                 push!(AB, latB2[i,:])
#             end
#             if distBB[1] < tol
#                 push!(BB, latB2[i,:])
#             end
#         end
#
#         latAA = transpose(hcat(AA...))
#         latBA = transpose(hcat(BA...))
#         latAB = transpose(hcat(AB...))
#         latBB = transpose(hcat(BB...))
#
#         try mkdir("data/"*ang_name*"_200k")
#         catch e
#         end
#
#         try write_lattice(latAA, "data/"*ang_name*"_200k/latticeAA.dat", max_radius)
#         catch e
#             println("AA: empty!")
#         end
#
#         try write_lattice(latBA, "data/"*ang_name*"_200k/latticeBA.dat", max_radius)
#         catch e
#             println("BA: empty!")
#         end
#
#         try write_lattice(latAB, "data/"*ang_name*"_200k/latticeAB.dat", max_radius)
#         catch e
#             println("AB: empty!")
#         end
#
#         try write_lattice(latBB, "data/"*ang_name*"_200k/latticeBB.dat", max_radius)
#         catch e
#             println("BB: empty!")
#         end
#
#         # WRITING POINTS OUT OF OVERLAP
#         try write_lattice(latA1, "data/"*ang_name*"_200k/latticeA1.dat", max_radius)
#         catch e
#             println("A1 lat is empty!")
#         end
#
#         try write_lattice(latB1, "data/"*ang_name*"_200k/latticeB1.dat", max_radius)
#         catch e
#             println("B1: empty!")
#         end
#
#         try write_lattice(latA2, "data/"*ang_name*"_200k/latticeA2.dat", max_radius)
#         catch e
#             println("A2 lat is empty!")
#         end
#
#         try write_lattice(latB2, "data/"*ang_name*"_200k/latticeB2.dat", max_radius)
#         catch e
#             println("B2 lat is empty!")
#         end
#
#         write_properties(p, q, j, steps, max_radius, a_top, a1_top, a2_top, a_bot, a1_bot, a2_bot, "data/"*ang_name*"_200k/properties.dat")
#     end
# end

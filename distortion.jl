include("hex_utils.jl")
using .HexUtils
using NearestNeighbors
using Printf


## finding smallest separation and corresponding angle
angles = Float64[]
AA_sep = Float64[]
BA_sep = Float64[]
AB_sep = Float64[]
BB_sep = Float64[]
AA_vec = []
BA_vec = []
AB_vec = []
BB_vec = []
AA_coord = []
BA_coord = []
AB_coord = []
BB_coord = []

open("AAstack_separations", "r") do f
    for line in eachline(f)
        if startswith(line, "0.")
            push!(angles, parse(Float64, line))
        elseif startswith(line, "AA:")
            aux = split(line, ":")
            push!(AA_sep, parse(Float64, aux[2]))
        elseif startswith(line, "BA:")
            aux = split(line, ":")
            push!(BA_sep, parse(Float64, aux[2]))
        elseif startswith(line, "AB:")
            aux = split(line, ":")
            push!(AB_sep, parse(Float64, aux[2]))
        elseif startswith(line, "BB:")
            aux = split(line, ":")
            push!(BB_sep, parse(Float64, aux[2]))
        elseif startswith(line, "AA_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AA_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BA_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BA_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "AB_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AB_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BB_vec:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BB_vec, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "AA_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AA_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BA_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BA_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "AB_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(AB_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        elseif startswith(line, "BB_coord:")
            aux = split(line, ":")
            aux2 = split(strip(aux[2]), ",")
            push!(BB_coord, [parse(Float64, lstrip(aux2[1], '[')), parse(Float64, aux2[2]), parse(Float64, rstrip(aux2[3], ']'))])
        end
    end
end

min_index = argmin(AB_sep)
println("Angle:       ", angles[min_index])
println("Separation:  ", AB_sep[min_index])
println("Sep. Vector: ", AB_vec[min_index])



## reading files path for the selected angle, loading lattice1
base_data_path = "data"
entries = readdir(base_data_path; join=true)
dataset_dirs = sort(filter(isdir, entries))
n_datasets = length(dataset_dirs)

path = dataset_dirs[min_index]
angle_name = path[findfirst("_",path)[1]+1:findlast("_", path)[1]-1]

latA1 = transpose(read_lattice_3d(path*"/latticeA1.dat"))
latB1 = transpose(read_lattice_3d(path*"/latticeB1.dat"))


## DISTORTION
## Method 2: changing the a1 and a2 vectors based on the separation vector
## and then generating the whole lattice again using modified lat vectors
## TODO: calculate area of primitive cell
angle, moire_period, max_radius, a1_top, a2_top, a1_bot, a2_bot = read_properties(path)

rotate_point!(a1_bot, angle, [0.0, 0.0])
rotate_point!(a2_bot, angle, [0.0, 0.0])

# cell_area_top = cell_area(a1_top, a2_top)
# cell_area_bot = cell_area(a1_bot, a2_bot)
# println(cell_area_top)
# println(cell_area_bot)

distance = sqrt(AB_coord[min_index][1]^2 + AB_coord[min_index][2]^2)
len_a = sqrt(a1_bot[1]^2 + a1_bot[2]^2)

cob_matrix = Matrix{Float64}(undef, 2, 2)
cob_matrix[1,1] = a1_top[1]
cob_matrix[1,2] = a2_top[1]
cob_matrix[2,1] = a1_top[2]
cob_matrix[2,2] = a2_top[2]
inv_cob = inv(cob_matrix)

coord_lat_basis = inv_cob * AB_coord[min_index][1:2]
# n_a1 = trunc(coord_lat_basis[1])
# n_a2 = trunc(coord_lat_basis[2])
n_a1 = coord_lat_basis[1]
n_a2 = coord_lat_basis[2]


println("Cardinal Coords: ", AB_coord[min_index])
println("Lat vecs Coords: ", [n_a1, n_a2])


println(a1_top)
println(a1_bot)




## finding AA, AB, BA, BB points for the angle and new A2, B2 lattices
# treeA1 = KDTree(latA1)
# treeB1 = KDTree(latB1)
#
# latA1 = transpose(latA1)
# latB1 = transpose(latB1)
# latA2 = transpose(latA2_distorted)
# latB2 = transpose(latB2_distorted)
#
# tol = 1.0e-3
# println("Tolerance:        ", tol)
# name = @sprintf("%6.4f", tol)
#
# AA = []
# BA = []
# AB = []
# BB = []
#
# n = minimum([size(latA1)[1], size(latB1)[1], size(latA2)[1], size(latB2)[1]])
#
# for i in 1:n
#     indAA, distAA = knn(treeA1, latA2[i,:], 1)
#     indBA, distBA = knn(treeB1, latA2[i,:], 1)
#     indAB, distAB = knn(treeA1, latB2[i,:], 1)
#     indBB, distBB = knn(treeB1, latB2[i,:], 1)
#     if distAA[1] < tol
#         push!(AA, latA2[i,:])
#     end
#     if distBA[1] < tol
#         push!(BA, latA2[i,:])
#     end
#     if distAB[1] < tol
#         push!(AB, latB2[i,:])
#     end
#     if distBB[1] < tol
#         push!(BB, latB2[i,:])
#     end
# end
#
# latAA = transpose(hcat(AA...))
# latBA = transpose(hcat(BA...))
# latAB = transpose(hcat(AB...))
# latBB = transpose(hcat(BB...))
#
#
# max_radius = maximum(latA1) - 10.0
# try write_lattice(latAA, path*"/latticeAA.dat", max_radius)
# catch e
#     println("AA lattice is empty!")
#     # println(e)
# end
#
# try write_lattice(latBA, path*"/latticeBA.dat", max_radius)
# catch e
#     println("BA lattice is empty!")
# end
#
# try write_lattice(latAB, path*"/latticeAB.dat", max_radius)
# catch e
#     println("AB lattice is empty!")
# end
#
# try write_lattice(latBB, path*"/latticeBB.dat", max_radius)
# catch e
#     println("BB lattice is empty!")
# end

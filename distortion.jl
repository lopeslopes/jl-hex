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
        end
    end
end

min_index = argmin(AB_sep)
println("Angle:      ", angles[min_index])
println("Separation: ", AB_sep[min_index])


## TODO: generate distorted lattices A2 and B2, based on the found separation



## finding AA, AB, BA, BB points for the angle and new A2, B2 lattices
base_data_path = "data"
entries = readdir(base_data_path; join=true)
dataset_dirs = sort(filter(isdir, entries))
n_datasets = length(dataset_dirs)

path = dataset_dirs[min_index]
angle_name = path[findfirst("_",path)[1]+1:findlast("_", path)[1]-1]
println(angle_name)

latA1 = transpose(read_lattice_3d(path*"/latticeA1.dat"))
latB1 = transpose(read_lattice_3d(path*"/latticeB1.dat"))
latA2 = transpose(read_lattice_3d(path*"/latticeA2.dat"))
latB2 = transpose(read_lattice_3d(path*"/latticeB2.dat"))

treeA1 = KDTree(latA1)
treeB1 = KDTree(latB1)

latA1 = transpose(latA1)
latB1 = transpose(latB1)
latA2 = transpose(latA2)
latB2 = transpose(latB2)

tol = 5.0e-3
println("Tolerance:        ", tol)
name = @sprintf("%6.4f", tol)

AA = []
BA = []
AB = []
BB = []

n = minimum([size(latA1)[1], size(latB1)[1], size(latA2)[1], size(latB2)[1]])

for i in 1:n
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


max_radius = maximum(latA1) - 10.0
try write_lattice(latAA, path*"/latticeAA.dat", max_radius)
catch e
    println("AA lattice is empty!")
    # println(e)
end

try write_lattice(latBA, path*"/latticeBA.dat", max_radius)
catch e
    println("BA lattice is empty!")
end

try write_lattice(latAB, path*"/latticeAB.dat", max_radius)
catch e
    println("AB lattice is empty!")
end

try write_lattice(latBB, path*"/latticeBB.dat", max_radius)
catch e
    println("BB lattice is empty!")
end

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

## reading files path for the selected angle, loading lattice1
# base_data_path = "data"
# entries = readdir(base_data_path; join=true)
# dataset_dirs = sort(filter(isdir, entries))
# n_datasets = length(dataset_dirs)
#
# path = dataset_dirs[min_index]

AB_stacking = false
angle = angles[min_index]
ang_name = @sprintf("%9.7f", angle)
if AB_stacking
    ang_name = "bernal_"*ang_name
else
    ang_name = "AA_"*ang_name
end
path = "data/"*ang_name*"_2M"

angle_name = path[findfirst("_",path)[1]+1:findlast("_", path)[1]-1]

## DISTORTION
## Method 2: changing the a1 and a2 vectors based on the separation vector
## and then generating the whole lattice again using modified lat vectors
angle, moire_period, max_radius, a1_top, a2_top, a1_bot, a2_bot = read_properties(path)

cob_matrix = Matrix{Float64}(undef, 2, 2)
cob_matrix[1,1] = a1_top[1]
cob_matrix[1,2] = a2_top[1]
cob_matrix[2,1] = a1_top[2]
cob_matrix[2,2] = a2_top[2]
inv_cob = inv(cob_matrix)

coord_lat_basis = inv_cob * AB_coord[min_index][1:2]
sep_vec_decomp = inv_cob * AB_vec[min_index][1:2]
m_coord = trunc(coord_lat_basis[1])
n_coord = trunc(coord_lat_basis[2])

s1 = sep_vec_decomp[1]/m_coord
s2 = sep_vec_decomp[2]/n_coord
sf = [s1, s2]

alpha1 = [1+s1, 0.0]
alpha2 = [0.0, 1+s2]

alpha1_cart = cob_matrix * alpha1
alpha2_cart = cob_matrix * alpha2

println("a1 cartesian: ", a1_top)
println("a2 cartesian: ", a2_top)

println("Alpha1 cartesian: ", alpha1_cart)
println("Alpha2 cartesian: ", alpha2_cart)

len_alpha1 = sqrt(alpha1_cart[1]^2 + alpha1_cart[2]^2)
len_alpha2 = sqrt(alpha2_cart[1]^2 + alpha2_cart[2]^2)
println("Alpha1 length: ", len_alpha1)
println("Alpha2 length: ", len_alpha2)

## MAKING NEW LATTICES BASED ON THE NEW ALPHA LATTICE VECTORS
n = 2000000
latA2 = read_lattice(path*"/latticeA2.dat")
latB2 = read_lattice(path*"/latticeB2.dat")
latA1_distorted = zeros(n ÷ 2, 2)
latB1_distorted = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1_distorted, latB1_distorted, alpha1_cart, alpha2_cart, false)

# cell_area_bot = cell_area(a1_bot, a2_bot)
# cell_area_top = cell_area(alpha1_cart, alpha2_cart)
# println("Cell area non-distorted: ", cell_area_bot)
# println("Cell area distorted:     ", cell_area_top)

# WRITING POINTS OUT OF OVERLAP
try write_lattice(latA1_distorted, path*"/latticeA1_dist.dat", max_radius)
catch e
    println(e)
end

try write_lattice(latB1_distorted, path*"/latticeB1_dist.dat", max_radius)
catch e
    println(e)
end

println(max_radius)

## TEST SECTION: READING LATTICE FROM FILE INSTEAD OF USING THE COMPUTED VALUES
## JUST TO BE EQUAL TO THE ACQUISITION OF LAT_A2 AND LAT_B2
latA1_distorted = 0
latB1_distorted = 0
# gc()
latA1_distorted = read_lattice(path*"/latticeA1_dist.dat", max_radius)
latB1_distorted = read_lattice(path*"/latticeB1_dist.dat", max_radius)


# finding AA, AB, BA, BB points for the angle and new A2, B2 lattices
treeA1 = KDTree(transpose(latA1_distorted))
treeB1 = KDTree(transpose(latB1_distorted))

# latA1_distorted = transpose(latA1_distorted)
# latB1_distorted = transpose(latB1_distorted)

tol = 1.0e-3
println("Tolerance:        ", tol)
name = @sprintf("%6.4f", tol)

AA = []
BA = []
AB = []
BB = []

n = minimum([size(latA1_distorted)[1], size(latB1_distorted)[1], size(latA2)[1], size(latB2)[1]])

println([size(latA1_distorted)[1], size(latB1_distorted)[1], size(latA2)[1], size(latB2)[1]])

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


max_radius = maximum(latA1_distorted) - 10.0

try write_lattice(latAA, path*"/latticeAA_dist.dat", max_radius)
catch e
    # println("AA lattice is empty!")
    println(e)
end

try write_lattice(latBA, path*"/latticeBA_dist.dat", max_radius)
catch e
    # println("BA lattice is empty!")
    println(e)
end

try write_lattice(latAB, path*"/latticeAB_dist.dat", max_radius)
catch e
    # println("AB lattice is empty!")
    println(e)
end

try write_lattice(latBB, path*"/latticeBB_dist.dat", max_radius)
catch e
    # println("BB lattice is empty!")
    println(e)
end

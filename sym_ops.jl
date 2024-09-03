include("hex_utils.jl")
using .HexUtils
using PyCall
using PyCall: LinearAlgebra
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Combinatorics
using StatsBase

ase = pyimport("ase")
spglib = pyimport("spglib")

lattice = [[4428.000519758899, 2484.072918780167, 0.0],
           [-4365.2705123956675, 2592.7244786925335, 0.0],
           [0.0, 0.0, 10000.0]]

a1_norm = LinearAlgebra.norm(lattice[1,:])
a2_norm = LinearAlgebra.norm(lattice[2,:])
a3_norm = LinearAlgebra.norm(lattice[3,:])

println(a1_norm)
println(a2_norm)
println(a3_norm)

positions = [[0.0, 0.0, 0.0]]
numbers = [1]

graphene = (lattice, positions, numbers)
sym_data = spglib.get_symmetry_dataset(graphene)
println("Group number:           ", sym_data.number)
println("Hall number:            ", sym_data.hall_number)
println("International notation: ", sym_data.international)
println("Hall notation:          ", sym_data.hall)
n_sym = size(sym_data.rotations)[1]
println("Number of symmetry operations of group: ", n_sym)

ax1 = subplot(111, projection="3d")

p1 = lattice[1][:]
p2 = lattice[2][:]
p3 = lattice[3][:]

ext_cell = hcat([p1,
                 p2,
                 -p1,
                 -p2,
                 [0.0, a1_norm, 0.0],
                 [0.0, -a1_norm, 0.0]]...)

ax1.scatter(ext_cell[1,:], ext_cell[2,:], ext_cell[3,:], s=200)
ext_cell = transpose(ext_cell)

aux_vec1 = p1/a1_norm
aux_vec2 = p2/a2_norm
aux_vec3 = ext_cell[5,:] + ext_cell[1,:]
aux_vec3 = aux_vec3/LinearAlgebra.norm(aux_vec3)
aux_vec4 = ext_cell[5,:] + ext_cell[2,:]
aux_vec4 = aux_vec4/LinearAlgebra.norm(aux_vec4)

grp_chr_names = ["" for i in 1:n_sym+1]
grp_chr = zeros(Int, n_sym+1)

for (i, rot) in enumerate(sym_data.rotations[:,3,3])
    # DISSECTING ROTATION MATRICES
    values, vectors = LinearAlgebra.eigen(rot)
    vectors = transpose(vectors)
    trc = LinearAlgebra.tr(rot)
    det = LinearAlgebra.det(rot)

    rot_axis = [0.0, 0.0, 1.0]
    rot_angle = pi/3.0

    gen_rot = zeros(Float64, 3, 3)
    for m in 1:3
        for l in 1:3
            aux_set = Set([1, 2, 3])
            delete!(aux_set, m)
            delete!(aux_set, l)
            n = first(aux_set)
            gen_rot[m, l] = (m == l) * cos(rot_angle) +
                            (det - cos(rot_angle)) * rot_axis[m] * rot_axis[l] -
                            sin(rot_angle) * levicivita([m, l, n]) * rot_axis[n]
        end
    end

    new_points = []
    for (jj, pn) in enumerate(ext_cell[:,3])
        # POINT
        new_p = gen_rot * pn
        push!(new_points, real(new_p))

        # SPIN TESTING
        spin_test = pn .+ [0.0, 0.0, 1.0]
        new_spin_test = real(gen_rot * spin_test)
        ds = new_spin_test .- new_p
        gc = 0
        if new_p in ext_cell
            gc = gc + Int(ds[3])
        end
        grp_chr[i] = grp_chr[i] + gc
    end

    new_points = hcat(new_points...)
    ax1.scatter(new_points[1, :], new_points[2, :], new_points[3, :])
end

println(grp_chr_names)

# count_dict = countmap(grp_chr_names)
# names_unique, counts = countmap(grp_chr_names)
# chr_unique = zeros(Int, length(names_unique))
#
# for (j, nu) in enumerate(names_unique)
#     ind = [i for i in 1:length(grp_chr_names) if grp_chr_names[i] == nu]
#     chr_final = 0
#     cont = 0
#     for i in ind
#         if grp_chr[i] != chr_final
#             chr_final = grp_chr[i]
#             cont += 1
#         end
#     end
#     chr_unique[j] = chr_final
# end
#
# println("----------------------------")
# println("CHARACTERES")
# println(rpad("sym_op", 7), rpad("chr", 3), rpad("mult", 4))
# for (n, c) in zip(names_unique, chr_unique)
#     mult = counts[findfirst(x -> x == n, names_unique)]
#     println(rpad(n, 7), rpad(string(c), 3), rpad(string(mult), 4))
# end

ax1.quiver(0.0, 0.0, 0.0, a1_norm*aux_vec1[1], a1_norm*aux_vec1[2], a1_norm*aux_vec1[3], color="red")
ax1.quiver(0.0, 0.0, 0.0, a1_norm*aux_vec2[1], a1_norm*aux_vec2[2], a1_norm*aux_vec2[3], color="red")
ax1.quiver(0.0, 0.0, 0.0, a1_norm*aux_vec3[1], a1_norm*aux_vec3[2], a1_norm*aux_vec3[3], color="red")
ax1.quiver(0.0, 0.0, 0.0, a1_norm*aux_vec4[1], a1_norm*aux_vec4[2], a1_norm*aux_vec4[3], color="red")

ax1.set_box_aspect((1, 1, 1))
ax1.set_xlim([-10000.0, 10000.0])
ax1.set_ylim([-10000.0, 10000.0])
ax1.set_zlim([-10000.0, 10000.0])

show()

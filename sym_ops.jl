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

spglib = pyimport("spglib")

lattice = [[4428.000519758899, 2484.072918780167, 0.0],
           [-4365.2705123956675, 2592.7244786925335, 0.0],
           [0.0, 0.0, 10000.0]]

# lattice = [ [sqrt(3)/2, 1/2, 0.0],
#             [-sqrt(3)/2, 1/2, 0.0],
#             [0.0, 0.0, 10000]]

a1_norm = LinearAlgebra.norm(lattice[1, :])
a2_norm = LinearAlgebra.norm(lattice[2, :])
a3_norm = LinearAlgebra.norm(lattice[3, :])

positions = [[0.0, 0.0, 0.0]]
numbers = [1]

graphene = (lattice, positions, numbers)
sym_data = spglib.get_symmetry_dataset(graphene)
n_sym = size(sym_data.rotations)[1]
println("Group number:           ", sym_data.number)
println("Hall number:            ", sym_data.hall_number)
println("International notation: ", sym_data.international)
println("Hall notation:          ", sym_data.hall)
println("Symmetry operations:    ", n_sym)

ax1 = subplot(111, projection="3d")

p1 = lattice[1][:]
rot60 = transpose(hcat([[1/2, -sqrt(3)/2, 0.0],
                        [sqrt(3)/2, 1/2, 0.0],
                        [0.0, 0.0, 1.0]]...))
p2 = rot60 * p1
p3 = rot60 * p2
p4 = rot60 * p3
p5 = rot60 * p4
p6 = rot60 * p5
ext_cell = hcat([p1,
                 p2,
                 p3,
                 p4,
                 p5,
                 p6]...)

ax1.scatter(ext_cell[1, :], ext_cell[2, :], ext_cell[3, :], s=200)
cell_tree = KDTree(ext_cell)
ext_cell = transpose(ext_cell)
n_pts = size(ext_cell)[1]

aux_vec1 = p1 / a1_norm
aux_vec2 = p3 / a1_norm
aux_vec3 = (p2 + p1) / LinearAlgebra.norm(p2 + p1)
aux_vec4 = (p2 + p3) / LinearAlgebra.norm(p2 + p3)
aux_vec5 = (p1 + p6) / LinearAlgebra.norm(p1 + p6)
aux_vec6 = p2 / a1_norm
aux_axis = [aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6]

grp_chr_names = ["" for i in 1:n_sym]
grp_chr = zeros(Int, n_sym)

for i in 1:n_sym
    rot = sym_data.rotations[i,:,:]
    det = LinearAlgebra.det(rot)

    rot_axis, rot_angle = analyze_sym_op!(rot, det, grp_chr_names, i, aux_axis)

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
    for jj in 1:n_pts
        # POINT
        pn = ext_cell[jj,1:3]
        new_p = gen_rot * pn
        push!(new_points, real(new_p))

        # SPIN TESTING
        spin_test = pn .+ [0.0, 0.0, 1.0]
        new_spin_test = gen_rot * spin_test
        ds = new_spin_test .- new_p
        gc = 0
        for pt_cell in eachrow(ext_cell)
            if (isapprox(new_p, pt_cell))
                gc = gc + Int(round(ds[3]))
            end
        end
        grp_chr[i] = grp_chr[i] + gc
    end

    new_points = hcat(new_points...)
    ax1.scatter(new_points[1, :], new_points[2, :], new_points[3, :])
    new_points = transpose(new_points)
end

count_dict = countmap(grp_chr_names)
chr_unique = zeros(Int, length(count_dict))
names_unique = ["" for i in 1:length(count_dict)]

for (j, op) in enumerate(keys(count_dict))
    ind = findall(name==op for name in grp_chr_names)
    names_unique[j] = op
    chr_final = 0
    cont = 0
    for i in ind
        if grp_chr[i] != chr_final
            chr_final = grp_chr[i]
            cont += 1
        end
    end
    chr_unique[j] = chr_final
end

println("----------------------------")
println("CHARACTERES")
println(lpad("sym_op", 9), lpad("chr", 5), lpad("mult", 6))
for (n, c) in zip(names_unique, chr_unique)
    mult = count_dict[n]
    println(lpad(n, 9), lpad(string(c), 5), lpad(string(mult), 6))
end

ax1.quiver(0.0, 0.0, 0.0, a1_norm * aux_vec1[1], a1_norm * aux_vec1[2], a1_norm * aux_vec1[3], color="red")
ax1.quiver(0.0, 0.0, 0.0, a1_norm * aux_vec2[1], a1_norm * aux_vec2[2], a1_norm * aux_vec2[3], color="red")
ax1.quiver(0.0, 0.0, 0.0, a1_norm * aux_vec3[1], a1_norm * aux_vec3[2], a1_norm * aux_vec3[3], color="red")
ax1.quiver(0.0, 0.0, 0.0, a1_norm * aux_vec4[1], a1_norm * aux_vec4[2], a1_norm * aux_vec4[3], color="red")

ax1.set_box_aspect((1, 1, 1))
ax1.set_xlim([-10000.0, 10000.0])
ax1.set_ylim([-10000.0, 10000.0])
ax1.set_zlim([-10000.0, 10000.0])

show()

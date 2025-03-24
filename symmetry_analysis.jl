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

latAA = read_lattice_3d("data/0.0192444/latticeAA.dat")
latBA = read_lattice_3d("data/0.0192444/latticeBA.dat")
latAB = read_lattice_3d("data/0.0192444/latticeAB.dat")
latBB = read_lattice_3d("data/0.0192444/latticeBB.dat")

origin = [0.0, 0.0, 0.0]
tree = KDTree(transpose(latAB))
ind_nbors, dist_nbors = knn(tree, origin, 7)
neighbors = [latAB[i,:] for i in ind_nbors]

a1 = neighbors[1][:]
aux = [a1[1], a1[2]]
origin_2d = [origin[1], origin[2]]
rotate_point!(aux, (2.0/3.0)*pi, origin_2d)
ind_a2, dist_a2 = nn(tree, [aux[1], aux[2], 0.0])
a2 = latAB[ind_a2,:]

lattice =  [[a1[1], a1[2], a1[3]],
            [a2[1], a2[2], a2[3]],
            [  0.0,   0.0, 1000000.0]]

positions = [[0.0, 0.0, 0.0]]
numbers = [1]

graphene = (lattice, positions, numbers)
sym_data = spglib.get_symmetry_dataset(graphene)

n_rot = size(sym_data.rotations)[1]
n_trans = size([sym_data.translations[i,:] for i in 1:size(sym_data.translations)[1] if sym_data.translations[i,:] != [0.0, 0.0, 0.0]])[1]

println("Group number:           ", sym_data.number)
println("Hall number:            ", sym_data.hall_number)
println("International notation: ", sym_data.international)
println("Rotations:              ", n_rot)
println("Translations:           ", n_trans)

grp_chr_names = ["" for i in 1:n_rot]
grp_chr = zeros(Int, n_rot)

# ax1 = subplot(111, projection="3d")

for i in 1:n_rot
    rot = sym_data.rotations[i,:,:]
    det = LinearAlgebra.det(rot)

    rot_axis, rot_angle = analyze_sym_op!(rot, grp_chr_names, i, lattice)

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
    op_viz = []
    for jj in 1:7
        # POINT
        pn = neighbors[jj][1:3]
        new_p = gen_rot * pn
        push!(new_points, real(new_p))
        
        # SPIN TESTING
        spin_test = pn .+ [0.0, 0.0, 1.0]
        new_spin_test = gen_rot * spin_test
        ds = new_spin_test .- new_p
        gc = 0
        for pt_cell in neighbors
            if (isapprox(new_p, pt_cell))
                gc = gc + Int(round(ds[3]))
            end
        end
        # if (isapprox(new_p, pn))
        #     gc = gc + 1 #Int(round(ds[3]))
        # end
        grp_chr[i] = grp_chr[i] + gc
    end

    new_points = hcat(new_points...)
    # ax1.scatter(new_points[1, :], new_points[2, :], new_points[3, :])
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

# ax1.set_box_aspect((1, 1, 1))
# ax1.set_xlim([-10.0, 10.0])
# ax1.set_ylim([-10.0, 10.0])
# ax1.set_zlim([-10.0, 10.0])
#
# show()

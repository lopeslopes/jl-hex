include("hex_utils.jl")
using .HexUtils
using PyCall
using PyCall: LinearAlgebra
pygui(:qt5)
using PyPlot
using Printf
using Combinatorics
using StatsBase
spglib = pyimport("spglib")


lattice = [[ 8742.840013282768, 5159.883256290869,     0.0],
           [-8840.009987150486, 4991.579924581863,     0.0],
           [               0.0,               0.0, 10000.0]]

p1 = [ 8742.840013282768,  5159.883256290869, 0.0]
p2 = [-97.16997386763752,  10151.46318086847, 0.0]
p3 = [-8840.009987150486,  4991.579924581863, 0.0]
p4 = [-8742.840013282766, -5159.883256290868, 0.0]
p5 = [  97.1699738677083, -10151.46318087216, 0.0]
p6 = [ 8840.009987150488, -4991.579924581928, 0.0]

a1_norm = LinearAlgebra.norm(lattice[1, :])
a2_norm = LinearAlgebra.norm(lattice[2, :])
a3_norm = LinearAlgebra.norm(lattice[3, :])

positions = [[0.0, 0.0, 0.0]]
numbers = [1]

graphene = (lattice, positions, numbers)
sym_data = spglib.get_symmetry_dataset(graphene, hall_number=0)
n_rot = size(sym_data.rotations)[1]
n_trans = size([sym_data.translations[i,:] for i in 1:size(sym_data.translations)[1] if sym_data.translations[i,:] != [0.0, 0.0, 0.0]])[1]

println("Group number:           ", sym_data.number)
println("Hall number:            ", sym_data.hall_number)
println("International notation: ", sym_data.international)
println("Hall notation:          ", sym_data.hall)
println("Rotations:              ", n_rot)
println("Translations:           ", n_trans)

ext_cell = hcat([p1, p2, p3, p4, p5, p6]...)
ax1 = subplot(111, projection="3d")
ax1.scatter(ext_cell[1, :], ext_cell[2, :], ext_cell[3, :], s=200)

ext_cell = transpose(ext_cell)
n_pts = size(ext_cell)[1]

aux_vec1 = p1 / a1_norm
aux_vec2 = p3 / a1_norm
aux_vec3 = (p2 + p1) / LinearAlgebra.norm(p2 + p1)
aux_vec4 = (p2 + p3) / LinearAlgebra.norm(p2 + p3)
aux_vec5 = (p1 + p6) / LinearAlgebra.norm(p1 + p6)
aux_vec6 = p2 / a1_norm
aux_axis = [aux_vec1, aux_vec2, aux_vec3, aux_vec4, aux_vec5, aux_vec6]

grp_chr_names = ["" for i in 1:n_rot]
grp_chr = zeros(Int, n_rot)

for i in 1:n_rot
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
    op_viz = []
    for jj in 1:n_pts
        # POINT
        pn = ext_cell[jj,1:3]
        new_p = gen_rot * pn
        push!(new_points, real(new_p))
        
        # VISUALIZATION
        max_viz = 50
        for k in 1:max_viz
            alpha_rot = rot_angle*real(k/max_viz)
            gen_rot2 = zeros(Float64, 3, 3)
            for m in 1:3
                for l in 1:3
                    aux_set = Set([1, 2, 3])
                    delete!(aux_set, m)
                    delete!(aux_set, l)
                    n = first(aux_set)
                    gen_rot2[m, l] = (m == l) * cos(alpha_rot) +
                                    (det - cos(alpha_rot)) * rot_axis[m] * rot_axis[l] -
                                    sin(alpha_rot) * levicivita([m, l, n]) * rot_axis[n]
                end
            end
            dgd_point = gen_rot2 * pn
            push!(op_viz, real(dgd_point))
        end
        op_viz = hcat(op_viz...)
        if (jj == 1)
            ax1.plot(op_viz[1,:], op_viz[2,:], op_viz[3,:])
        end
        op_viz = []

        # SPIN TESTING
        spin_test = pn .+ [0.0, 0.0, 1.0]
        new_spin_test = gen_rot * spin_test
        ds = new_spin_test .- new_p
        gc = 0
        if (isapprox(new_p, pn))
            gc = gc + 1 #Int(round(ds[3]))
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

ax1.set_box_aspect((1, 1, 1))
ax1.set_xlim([-10000.0, 10000.0])
ax1.set_ylim([-10000.0, 10000.0])
ax1.set_zlim([-10000.0, 10000.0])

show()

include("hex_utils.jl")
using .HexUtils
using PyCall
using PyCall: LinearAlgebra
pygui(:qt5)
using PyPlot
using Printf
using Combinatorics
using NearestNeighbors
using StatsBase
spglib = pyimport("spglib")


lattice = [[ 8742.840013282768, 5159.883256290869,     0.0],
           [-8840.009987150486, 4991.579924581863,     0.0],
           [               0.0,               0.0, 10000.0]]

norm_a1 = LinearAlgebra.norm(lattice[1, :])

positions = [[0.0, 0.0, 0.0]]
numbers = [1]
graphene = (lattice, positions, numbers)
sym_data = spglib.get_symmetry_dataset(graphene, hall_number=0)
n_rot = size(sym_data.rotations)[1]

grp_chr_names = ["" for i in 1:n_rot]
op_name = []
op_cartesian = []
op_latbase = []
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

    if (det == 1) || (det == -1)
        push!(op_cartesian, gen_rot)
        push!(op_latbase, rot)
        push!(op_name, grp_chr_names[i])
    end
end

latAA = transpose(read_lattice_3d("data/0.0191435_aa_stack/latticeAA.dat", norm_a1+200.0))
latBA = transpose(read_lattice_3d("data/0.0191435_aa_stack/latticeBA.dat", norm_a1+200.0))
latAB = transpose(read_lattice_3d("data/0.0191435_aa_stack/latticeAB.dat", norm_a1+200.0))
latBB = transpose(read_lattice_3d("data/0.0191435_aa_stack/latticeBB.dat", norm_a1+200.0))

tol = 1e-8

current_lat = latAA
tree = KDTree(current_lat)

ax1 = subplot(111)

num_out = 0
total_atoms = Int(size(current_lat)[2])
# println(lpad("ind", 3), "|", lpad("op_name", 8), "|", lpad("Atoms", 6), "|", lpad("Displaced", 10))
# for op_ind in 8:8
for op_ind in 1:size(op_cartesian)[1]
    operation = op_cartesian[op_ind][:,:]

    out_pt = []
    global num_out = 0
    for point in eachcol(current_lat)
        aux_p = operation * point
        ind_pt, dist_pt = nn(tree, aux_p)
        if (dist_pt > tol)
            global num_out += 1
            push!(out_pt, aux_p)
        end
    end
    # println(lpad(string(op_ind), 3), "|", lpad(op_name[op_ind], 8), "|", lpad(string(total_atoms), 6), "|", lpad(string(num_out), 10))
    if (num_out == 0)
        println(op_name[op_ind])
    end
    out_pt = hcat(out_pt...)
    try
        ax1.scatter(out_pt[1,:], out_pt[2,:], s=20, color=("red", 0.7))
    catch e
    end
end

# ax1.scatter(current_lat[1,:], current_lat[2,:], s=50, color=("green", 0.3))
#
# ax1.set_xlim([-10000, 10000])
# ax1.set_ylim([-10000, 10000])
# ax1.set_aspect("equal")
#
# show()

### AA STACKED:
# AA: todos
# BA: E, C3, S3, sigma_h
# AB: E, C3, S3, sigma_h
# BB: E, C3, S3, sigma_h, C''2, sigma_v
#
### BERNAL STACKED:
# AA: E, C3, S3, sigma_h
# BA: E, C3, S3, sigma_h, C'2, sigma_d
# AB: todos
# BB: E, C3, S3, sigma_h

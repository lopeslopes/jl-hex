module HexUtils

using LinearAlgebra

export create_honeycomb_lattice!, write_lattice, rotate_lattice!, read_lattice, read_lattice_3d, analyze_sym_op!


function create_honeycomb_lattice!(latticeA::Array{Float64,2}, latticeB::Array{Float64,2}, a::Float64, ab_stacking::Bool)
    num_columns = 2*div(isqrt(size(latticeA,1)*2),3)
    angle = (60.0/180.0) * π

    v1 = [a, 0.0]
    v2 = [a*cos(angle), a*sin(angle)]
    v3 = [-a*cos(angle), a*sin(angle)]

    d = sqrt((a^2)/(2.0*(1.0-cos(2.0*angle))))
    d1 = [d*cos(angle/2.0), d*sin(angle/2.0)]

    origin_a = [0.0, 0.0]
    origin_b = origin_a + d1

    row = 1
    i = 1
    while(i < size(latticeA,1))
        for j=1:num_columns
            if (i > size(latticeA,1))
                break
            end
            latticeA[i,:] = origin_a + Float64(j-1)*v1
            latticeB[i,:] = origin_b + Float64(j-1)*v1
            i = i + 1
        end
        row  = row + 1
        if (row % 2 == 1)
            origin_a = origin_a + v2
            origin_b = origin_b + v2
        else
            origin_a = origin_a + v3
            origin_b = origin_b + v3
        end
    end

    i0 = div((row ÷ 2) * num_columns + div(num_columns,2), 1)
    lat_origin = latticeA[i0,:]

    latticeA[:,1] = latticeA[:,1] .- lat_origin[1]
    latticeA[:,2] = latticeA[:,2] .- lat_origin[2]
    latticeB[:,1] = latticeB[:,1] .- lat_origin[1]
    latticeB[:,2] = latticeB[:,2] .- lat_origin[2]
    
    if (ab_stacking)
        latticeA[:,1] = latticeA[:,1] .- d1[1]
        latticeA[:,2] = latticeA[:,2] .- d1[2]
        latticeB[:,1] = latticeB[:,1] .- d1[1]
        latticeB[:,2] = latticeB[:,2] .- d1[2]
    end
end

function write_lattice(lattice, filename)
    open(filename, "w") do file
        for i = 1:size(lattice, 1)
            println(file, lattice[i, 1], ";", lattice[i, 2])
        end
    end
end

function read_lattice(filename, max_dist=0.0)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2])]
                dist = sqrt(aux_v[1]^2 + aux_v[2]^2)
                if ((max_dist == 0.0) || (dist < max_dist))
                    push!(lat, aux_v)
                end
            end
        end
    end
    lat = transpose(hcat(lat...))
    return lat
end

function read_lattice_3d(filename, max_dist=0.0)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2]), 0.0]
                dist = sqrt(aux_v[1]^2 + aux_v[2]^2 + aux_v[3]^2)
                if ((max_dist == 0.0) || (dist < max_dist))
                    push!(lat, aux_v)
                end
            end
        end
    end
    lat = transpose(hcat(lat...))
    return lat
end

function rotate_lattice!(lattice, angle, pivot)
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    for i = 1:size(lattice, 1)
        aux1 = lattice[i, :] - pivot
        aux1 = rot_matrix * aux1
        lattice[i, :] .= aux1 .+ pivot
    end
end

function analyze_sym_op!(rot_matrix, grp_chr_names, i, lattice_vecs)
    trc = LinearAlgebra.tr(rot_matrix)
    det = LinearAlgebra.det(rot_matrix)
    solution = LinearAlgebra.eigen(rot_matrix)
    values = solution.values
    vectors = transpose(solution.vectors)

    a1_norm = LinearAlgebra.norm(lattice_vecs[1, :])

    p_aux = (lattice_vecs[2] + lattice_vecs[1]) / LinearAlgebra.norm(lattice_vecs[2] + lattice_vecs[1])
    aux_v1 = lattice_vecs[1] / a1_norm
    aux_v2 = lattice_vecs[2] / a1_norm
    aux_v3 = (aux_v1 + p_aux) / LinearAlgebra.norm(aux_v1 + p_aux)
    aux_v4 = (aux_v2 + p_aux) / LinearAlgebra.norm(aux_v2 + p_aux)
    aux_v5 = (aux_v1 - aux_v2) / LinearAlgebra.norm(aux_v2 - aux_v1)
    aux_v6 = p_aux / LinearAlgebra.norm(p_aux)
    aux_axis = [aux_v1, aux_v2, aux_v3, aux_v4, aux_v5, aux_v6]

    ax_ind = findall(isapprox(v, det) for v in values)

    axis = real(vectors[ax_ind[end],1:3])
    angle = acos((trc - det)/2)

    if (det == 1.0)
        if (trc == 3)
            grp_chr_names[i] = "E"
        elseif (trc == 2)
            if ("C6" in grp_chr_names)
                angle = angle + 2*(pi - angle)
            end
            grp_chr_names[i] = "C6"
        elseif (trc == 0)
            if ("C3" in grp_chr_names)
                angle = angle + 2*(pi - angle)
            end
            grp_chr_names[i] = "C3"
        elseif (trc == -1)
            if (isapprox(axis, [1.0, 0.0, 0.0]))
                grp_chr_names[i] = "C*2"
                axis = aux_axis[5][:]
            elseif (isapprox(axis, [0.0, 1.0, 0.0]))
                grp_chr_names[i] = "C**2"
                axis = aux_axis[6][:]
            elseif (isapprox(axis, [0.0, 0.0, 1.0]))
                grp_chr_names[i] = "C2"
            elseif (isapprox(axis[1], axis[2]))
                axis = aux_axis[3][:]
                grp_chr_names[i] = "C*2"
            elseif (isapprox(axis[1], -axis[2]))
                axis = aux_axis[4][:]
                grp_chr_names[i] = "C*2"
            elseif (isapprox(axis[1], 2*axis[2]))
                axis = aux_axis[1][:]
                grp_chr_names[i] = "C**2"
            elseif (isapprox(2*axis[1], axis[2]))
                axis = aux_axis[2][:]
                grp_chr_names[i] = "C**2"
            end
        end
    elseif (det == -1.0)
        if (trc == -3)
            grp_chr_names[i] = "i"
        elseif (trc == -2)
            if ("S3" in grp_chr_names)
                angle = angle + 2*(pi - angle)
            end
            grp_chr_names[i] = "S3"
        elseif (trc == 0)
            if ("S6" in grp_chr_names)
                angle = angle + 2*(pi - angle)
            end
            grp_chr_names[i] = "S6"
        elseif (trc == 1)
            if (isapprox(axis, [1.0, 0.0, 0.0]))
                grp_chr_names[i] = "sigma_v"
                axis = aux_axis[5][:]
            elseif (isapprox(axis, [0.0, 1.0, 0.0]))
                grp_chr_names[i] = "sigma_d"
                axis = aux_axis[6][:]
            elseif (isapprox(axis, [0.0, 0.0, 1.0]))
                grp_chr_names[i] = "sigma_h"
            elseif (isapprox(axis[1], axis[2]))
                axis = aux_axis[3][:]
                grp_chr_names[i] = "sigma_v"
            elseif (isapprox(axis[1], -axis[2]))
                axis = aux_axis[4][:]
                grp_chr_names[i] = "sigma_v"
            elseif (isapprox(axis[1], 2*axis[2]))
                axis = aux_axis[1][:]
                grp_chr_names[i] = "sigma_d"
            elseif (isapprox(2*axis[1], axis[2]))
                axis = aux_axis[2][:]
                grp_chr_names[i] = "sigma_d"
            end
        end
    end

    if (grp_chr_names[i] == "")
        println(axis)
    end

    return axis, angle
end

end

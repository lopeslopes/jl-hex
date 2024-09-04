module HexUtils

using LinearAlgebra

export create_honeycomb_lattice!, create_honeycomb_lattice_fractional!, write_lattice, rotate_lattice!, read_lattice, rotate_3d, sym_op_d6h, analyze_sym_op!


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

function create_honeycomb_lattice_fractional!(latticeA::Array{Float64,2}, latticeB::Array{Float64,2}, a::Float64, ab_stacking::Bool)
    num_columns = 2*div(isqrt(size(latticeA,1)*2),3)
    #d = sqrt((a^2)/(2.0*(1.0-cos(2.0*angle))))
    #d1 = [d*cos(angle/2.0), d*sin(angle/2.0)]

    origin_a = [0.0, 0.0]
    origin_b = [0.0, 0.0] #origin_a + d1

    row = 1
    i = 1
    while(i < size(latticeA,1))
        for j=1:num_columns
            if (i > size(latticeA,1))
                break
            end
            latticeA[i,:] = [origin_a[1] + Float64(j-1), origin_a[2]]
            latticeB[i,:] = [origin_b[1] + Float64(j-1), origin_b[2]]
            i = i + 1
        end
        row  = row + 1
        if (row % 2 == 1)
            origin_a = origin_a + [0.0, 1.0]
            origin_b = origin_b + [0.0, 1.0]
        else
            origin_a = origin_a + [-1.0, 1.0]
            origin_b = origin_b + [-1.0, 1.0]
        end
    end

    i0 = div((row ÷ 2) * num_columns + div(num_columns,2), 1)
    lat_origin = latticeA[i0,:]

    #latticeA[:,1] = latticeA[:,1] .- lat_origin[1]
    #latticeA[:,2] = latticeA[:,2] .- lat_origin[2]
    #latticeB[:,1] = latticeB[:,1] .- lat_origin[1]
    #latticeB[:,2] = latticeB[:,2] .- lat_origin[2]
    
    #if (ab_stacking)
    #    latticeA[:,1] = latticeA[:,1] .- d1[1]
    #    latticeA[:,2] = latticeA[:,2] .- d1[2]
    #    latticeB[:,1] = latticeB[:,1] .- d1[1]
    #    latticeB[:,2] = latticeB[:,2] .- d1[2]
    #end
end

function write_lattice(lattice, filename)
    open(filename, "w") do file
        for i = 1:size(lattice, 1)
            println(file, lattice[i, 1], ";", lattice[i, 2])
        end
    end
end

function read_lattice(filename)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2])]
                push!(lat, aux_v)
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

function rotate_3d(point, angle, rot_axis)
    # ROTATION AXIS: 1=x, 2=y, 3=z
    point_3d = [point[1], point[2], 0.0]

    if rot_axis == 1
        rot_matrix = [1     0           0     ; 
                      0 cos(angle) -sin(angle); 
                      0 sin(angle)  cos(angle)]
    elseif rot_axis == 2
        rot_matrix = [ cos(angle) 0 sin(angle); 
                           0      1     0     ; 
                      -sin(angle) 0 cos(angle)]
    else
        rot_matrix = [cos(angle) -sin(angle) 0; 
                      sin(angle)  cos(angle) 0; 
                          0           0      1]
    end

    point_3d = rot_matrix * point_3d
    new_point = [point_3d[1], point_3d[2]]
    return new_point
end

function sym_op_d6h(point, operation)
    point_3d = [point[1], point[2], 0.0]

    if operation == 1 # E
        rot_matrix = [1 0 0; 
                      0 1 0; 
                      0 0 1]
    elseif operation == 2 # C6_1
        angle = pi/3
        rot_matrix = [cos(angle) -sin(angle) 0; 
                      sin(angle)  cos(angle) 0;
                          0           0      1]
    elseif operation == 3 # C6_2
        angle = 5*pi/3
        rot_matrix = [cos(angle) -sin(angle) 0; 
                      sin(angle)  cos(angle) 0;
                          0           0      1]
    elseif operation == 4 # C3_1
        angle = 2*pi/3
        rot_matrix = [cos(angle) -sin(angle) 0; 
                      sin(angle)  cos(angle) 0;
                          0           0      1]
    elseif operation == 5 # C3_2
        angle = 4*pi/3
        rot_matrix = [cos(angle) -sin(angle) 0; 
                      sin(angle)  cos(angle) 0;
                          0           0      1]
    elseif operation == 6 # C2
        rot_matrix = [-1  0  0; 
                       0 -1  0; 
                       0  0  1]
    elseif operation == 7 # C'2_1
        rot_matrix = [1  0  0; 
                      0 -1  0; 
                      0  0 -1]
    elseif operation == 8 # C'2_2
        angle = 4*pi/3
        rot_matrix = [ cos(angle) -sin(angle)  0; 
                      -sin(angle) -cos(angle)  0;
                           0           0      -1]
    elseif operation == 9 # C'2_3
        angle = 4*pi/3
        rot_matrix = [cos(angle)  sin(angle)  0; 
                      sin(angle) -cos(angle)  0;
                          0           0      -1]
    elseif operation == 10 # C''2_1
        angle = 4*pi/3
        rot_matrix = [-cos(angle) -sin(angle)  0; 
                      -sin(angle)  cos(angle)  0;
                           0           0      -1]
    elseif operation == 11 # C''2_2
        rot_matrix = [-1  0  0; 
                       0  1  0; 
                       0  0 -1]
    elseif operation == 12 # C''2_3
        angle = 4*pi/3
        rot_matrix = [-cos(angle) sin(angle)  0; 
                       sin(angle) cos(angle)  0;
                           0          0      -1]
    elseif operation == 13 # i
        rot_matrix = [-1  0  0; 
                       0 -1  0; 
                       0  0 -1]
    elseif operation == 14 # S3_1
        angle = 4*pi/3
        rot_matrix = [ cos(angle) sin(angle)  0; 
                      -sin(angle) cos(angle)  0;
                           0          0      -1]
    elseif operation == 15 # S3_2
        angle = 4*pi/3
        rot_matrix = [cos(angle) -sin(angle)  0; 
                      sin(angle)  cos(angle)  0;
                          0           0      -1]
    elseif operation == 16 # S6_1
        angle = 4*pi/3
        rot_matrix = [-cos(angle)  sin(angle)  0; 
                      -sin(angle) -cos(angle)  0;
                           0           0      -1]
    elseif operation == 17 # S6_2
        angle = 4*pi/3
        rot_matrix = [-cos(angle) -sin(angle)  0; 
                       sin(angle) -cos(angle)  0;
                           0           0      -1]
    elseif operation == 18 # sigma_h
        rot_matrix = [1  0  0; 
                      0  1  0; 
                      0  0 -1]
    elseif operation == 19 # sigma_d_1
        angle = 4*pi/3
        rot_matrix = [-cos(angle) -sin(angle)  0; 
                      -sin(angle)  cos(angle)  0;
                           0           0       1]
    elseif operation == 20 # sigma_d_2
        rot_matrix = [-1  0  0; 
                       0  1  0; 
                       0  0  1]
    elseif operation == 21 # sigma_d_3
        angle = 4*pi/3
        rot_matrix = [-cos(angle)  sin(angle)  0; 
                       sin(angle)  cos(angle)  0;
                           0           0       1]
    elseif operation == 22 # sigma_v_1
        rot_matrix = [1  0  0; 
                      0 -1  0; 
                      0  0  1]
    elseif operation == 23 # sigma_v_2
        angle = 4*pi/3
        rot_matrix = [ cos(angle) -sin(angle)  0; 
                      -sin(angle) -cos(angle)  0;
                           0           0       1]
    elseif operation == 24 # sigma_v_3
        angle = 4*pi/3
        rot_matrix = [cos(angle)  sin(angle)  0; 
                      sin(angle) -cos(angle)  0;
                          0           0       1]
    end

    point_3d = rot_matrix * point_3d
    new_point = [point_3d[1], point_3d[2]]
    return new_point
end

function analyze_sym_op!(rot_matrix, det, grp_chr_names, i, aux_axis)
    trc = LinearAlgebra.tr(rot_matrix)
    solution = LinearAlgebra.eigen(rot_matrix)
    values = solution.values
    vectors = transpose(solution.vectors)

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
            elseif (isapprox(axis, [0.0, 1.0, 0.0]))
                grp_chr_names[i] = "C**2"
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
            elseif (isapprox(axis, [0.0, 1.0, 0.0]))
                grp_chr_names[i] = "sigma_d"
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

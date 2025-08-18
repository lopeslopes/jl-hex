module HexUtils

using LinearAlgebra

export create_honeycomb_lattice!, write_lattice, rotate_lattice!, rotate_point!, read_lattice, read_lattice_3d, magic_angle, write_properties, read_properties, cell_area


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

function write_lattice(lattice, filename, max_dist=0.0, min_dist=0.0)
    open(filename, "w") do file
        for i = 1:size(lattice, 1)
            dist = sqrt(lattice[i, 1]^2 + lattice[i, 2]^2)
            if ((max_dist == 0.0) || (dist < max_dist))
                if ((min_dist == 0.0) || (dist > min_dist))
                    println(file, lattice[i, 1], ";", lattice[i, 2])
                end
            end
        end
    end
end

function read_lattice(filename, max_dist=0.0, min_dist=0.0)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2])]
                dist = sqrt(aux_v[1]^2 + aux_v[2]^2)
                if ((max_dist == 0.0) || (dist < max_dist))
                    if ((min_dist == 0.0) || (dist > min_dist))
                        push!(lat, aux_v)
                    end
                end
            end
        end
    end
    lat = transpose(hcat(lat...))
    return lat
end

function read_lattice_3d(filename, max_dist=0.0, min_dist=0.0)
    lat = []
    open(filename, "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, ";")
                aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2]), 0.0]
                dist = sqrt(aux_v[1]^2 + aux_v[2]^2 + aux_v[3]^2)
                if ((max_dist == 0.0) || (dist < max_dist))
                    if ((min_dist == 0.0) || (dist > min_dist))
                        push!(lat, aux_v)
                    end
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

function rotate_point!(point, angle, pivot)
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    aux1 = point .- pivot
    aux1 = rot_matrix * aux1
    point .= aux1 .+ pivot
end

function magic_angle(p, q)
    angle = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2))) 
end

function write_properties(p, q, i, steps, max_radius, filename)
    open(filename, "w") do file
        println(file, "p=", p)
        println(file, "q=", q)
        println(file, "i=", i)
        println(file, "steps=", steps)
        println(file, "max_radius=", max_radius)
    end
end

function read_properties(path)
    p = 0
    q = 0
    i = 0
    steps = 0
    max_radius = 0.0
    open(path*"/properties.dat", "r") do file
        data = readlines(file)
        for line in data
            if line != "\n"
                aux = split(line, "=")
                if aux[1] == "p"
                    aux_p = parse(Float64, aux[2])
                    p = trunc(Int, aux_p)
                elseif aux[1] == "q"
                    aux_q = parse(Float64, aux[2])
                    q = trunc(Int, aux_q)
                elseif aux[1] == "i"
                    aux_i = parse(Float64, aux[2])
                    i = trunc(Int, aux_i)
                elseif aux[1] == "steps"
                    aux_steps = parse(Float64, aux[2])
                    steps = trunc(Int, aux_steps)
                elseif aux[1] == "max_radius"
                    aux_radius = parse(Float64, aux[2])
                    max_radius = aux_radius
                end
            end
        end
    end
    angle_i = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
    angle_f = acos((3.0*((q-1)^2) - (p^2))/(3.0*((q-1)^2) + (p^2)))
    angle = angle_i + (i/steps)*(angle_f  - angle_i)

    a = 2.46
    moire_period = a/(2*sin(angle/2))

    return angle, moire_period, max_radius
end

function cell_area(a1, a2)
    a1_3d = [a1[1], a1[2], 0.0]
    a2_3d = [a2[1], a2[2], 0.0]
    area_lat_cell = LinearAlgebra.norm(LinearAlgebra.cross(a1_3d, a2_3d))
    return area_lat_cell
end

end

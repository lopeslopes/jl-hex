module HexUtils
export create_honeycomb_lattice!, write_lattice, rotate_lattice!

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
        latticeA[:,1] = latticeA[:,1] .+ d1[1]
        latticeA[:,2] = latticeA[:,2] .+ d1[2]
        latticeB[:,1] = latticeB[:,1] .+ d1[1]
        latticeB[:,2] = latticeB[:,2] .+ d1[2]
    end
end

function write_lattice(lattice, filename)
    open(filename, "w") do file
        for i = 1:size(lattice, 1)
            println(file, lattice[i, 1], ";", lattice[i, 2])
        end
    end
end

function rotate_lattice!(lattice, angle, pivot)
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)]

    for i = 1:size(lattice, 1)
        aux1 = lattice[i, :] - pivot
        aux1 = rot_matrix * aux1
        lattice[i, :] .= aux1 .+ pivot
    end
end

end

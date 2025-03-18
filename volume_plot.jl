using PyCall
pygui(:qt5)
using PyPlot


angle = []
vol_ws = []
vol_pv = []
open("cell_volume.dat", "r") do file
    data = readlines(file)
    for line in data
        if (line != data[1]) && (line != "\n") && (line != "")
            aux = split(line, ";")
            aux_v = [parse(Float64, aux[1]), parse(Float64, aux[2]), parse(Float64, aux[3])]
            push!(angle, aux_v[1])
            push!(vol_ws, aux_v[2])
            push!(vol_pv, aux_v[3])
        end
    end
end

ax1 = subplot(111)
ax1.scatter(angle, vol_ws)
ax1.scatter(angle, vol_pv)

show()

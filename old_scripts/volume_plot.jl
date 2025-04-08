using WGLMakie


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


with_theme(theme_latexfonts()) do
    fig = Figure()
    ax1 = Axis(fig[1, 1],
               title = "Area of unit cell",
               xlabel = L"\theta\, (\text{rad})",
               ylabel = L"\text{Area}\, (\AA^2)")

    scatter!(ax1, angle, vol_pv, label="Volume")
    # axislegend(position = :rb)
    save("test_plot.png", fig)
end

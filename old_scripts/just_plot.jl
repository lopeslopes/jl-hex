include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using Printf


#tol = 5.0e-3
tol = 1.0e-2

name1 = @sprintf("%6.4f", 1e-2)
name5 = @sprintf("%6.4f", 5e-3)

dirlist = readdir("data/")
#println(dirlist)

list_001 = []
list_0005 = []
list_bernal_1 = []
list_bernal_5 = []
for folder in dirlist
    if sizeof(folder) == 9
        push!(list_0005, folder)
    elseif sizeof(folder) == 14
        push!(list_001, folder)
    elseif sizeof(folder) == 21
        push!(list_bernal_1, folder)
    elseif sizeof(folder) == 16
        push!(list_bernal_5, folder)
    end
end

#println(list_bernal_1, list_bernal_5)

for (folder, folder5) in zip(list_001, list_0005)
    latAA = read_lattice("data/"*folder*"/200M_"*name1*"_AA.dat")
    latBA = read_lattice("data/"*folder*"/200M_"*name1*"_BA.dat")
    latAB = read_lattice("data/"*folder*"/200M_"*name1*"_AB.dat")
    latBB = read_lattice("data/"*folder*"/200M_"*name1*"_BB.dat")
    
    latAA5 = read_lattice("data/"*folder5*"/200M_"*name5*"_AA.dat")
    latBA5 = read_lattice("data/"*folder5*"/200M_"*name5*"_BA.dat")
    latAB5 = read_lattice("data/"*folder5*"/200M_"*name5*"_AB.dat")
    latBB5 = read_lattice("data/"*folder5*"/200M_"*name5*"_BB.dat")


    ax1 = subplot(121,aspect=1)
    try ax1.scatter(latAA[:,1], latAA[:,2], s=2, color="blue")
    catch e
        println("AA lattice is empty!")
    end
    try ax1.scatter(latBA[:,1], latBA[:,2], s=2, color="green")
    catch e
        println("BA lattice is empty!")
    end
    try ax1.scatter(latAB[:,1], latAB[:,2], s=2, color="orange")
    catch e
        println("AB lattice is empty!")
    end
    try ax1.scatter(latBB[:,1], latBB[:,2], s=2, color="red")
    catch e
        println("BB lattice is empty!")
    end
    
    ax1.set_xlim([-9000, 9000])
    ax1.set_ylim([-9000, 9000])

    #ax1.yaxis.set_label_position("right")
    #ax1.yaxis.tick_right()
    #ax1.tick_params(labelbottom=false, bottom=false)
    ax1.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)
    #legend(["AA", "BA", "AB", "BB"])


    ax2 = subplot(122,aspect=1)
    try ax2.scatter(latAA5[:,1], latAA5[:,2], s=2, color="blue")
    catch e
        println("AA lattice is empty!")
    end
    try ax2.scatter(latBA5[:,1], latBA5[:,2], s=2, color="green")
    catch e
        println("BA lattice is empty!")
    end
    try ax2.scatter(latAB5[:,1], latAB5[:,2], s=2, color="orange")
    catch e
        println("AB lattice is empty!")
    end
    try ax2.scatter(latBB5[:,1], latBB5[:,2], s=2, color="red")
    catch e
        println("BB lattice is empty!")
    end
    
    ax2.set_xlim([-9000, 9000])
    ax2.set_ylim([-9000, 9000])

    ax2.yaxis.set_label_position("right")
    ax2.yaxis.tick_right()
    ax2.tick_params(labelbottom=false, bottom=false)
    ax1.set_ylabel(L"$\AA$")
    #ax2.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)
    #legend(["AA", "BA", "AB", "BB"])
    
    title1 = folder[1:9]
    angle1 = parse(Float64, title1)
    ang_deg = angle1*180/pi
    title2 = @sprintf("%6.4f", ang_deg)
    ax1.set_title(L"$\theta=$"*title2*L"$^{\circ}$, $\delta=0.01 \AA$") 
    ax2.set_title(L"$\theta=$"*title2*L"$^{\circ}$, $\delta=0.005 \AA$")

    legend(["AA", "BA", "AB", "BB"], loc="center right", bbox_to_anchor=(0.1, 0.5), fancybox=true, shadow=true, markerscale=5)

    savefig("results/angvar_AA/"*folder5*".png", format="png", dpi=300, bbox_inches="tight")
    clf()
end
#show()

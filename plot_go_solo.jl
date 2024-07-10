include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using Printf

latAA1 = read_lattice("data/0.0191434/0.0050_AA.dat")
latBA1 = read_lattice("data/0.0191434/0.0050_BA.dat")
latAB1 = read_lattice("data/0.0191434/0.0050_AB.dat")
latBB1 = read_lattice("data/0.0191434/0.0050_BB.dat")

ax1 = subplot(111,aspect=1)
try ax1.scatter(latAA1[:,1], latAA1[:,2], s=10, color="c")
catch e
    println("AA lattice is empty!")
end
try ax1.scatter(latBA1[:,1], latBA1[:,2], s=10, color="green")
catch e
    println("BA lattice is empty!")
end
try ax1.scatter(latAB1[:,1], latAB1[:,2], s=10, color="orange")
catch e
    println("AB lattice is empty!")
end
try ax1.scatter(latBB1[:,1], latBB1[:,2], s=10, color="red")
catch e
    println("BB lattice is empty!")
end
ax1.set_xlim([-3500, 3500])
ax1.set_ylim([-3500, 3500])
#legend(["AA", "BA", "AB", "BB"])
ax1.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)

ax1.set_title(L"$\theta=1.0968^{\circ}$, $\delta=0.005 \AA$") 

legend(["AA", "BA", "AB", "BB"], loc="center left", bbox_to_anchor=(1, 0.5))

show()

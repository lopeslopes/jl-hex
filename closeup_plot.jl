include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using Printf

latAA1 = read_lattice("data/0.0189290_0.01/200M_0.0100_AA.dat")
latBA1 = read_lattice("data/0.0189290_0.01/200M_0.0100_BA.dat")
latAB1 = read_lattice("data/0.0189290_0.01/200M_0.0100_AB.dat")
latBB1 = read_lattice("data/0.0189290_0.01/200M_0.0100_BB.dat")

latAA2 = read_lattice("data/0.0195706_0.01/200M_0.0100_AA.dat")
latBA2 = read_lattice("data/0.0195706_0.01/200M_0.0100_BA.dat")
latAB2 = read_lattice("data/0.0195706_0.01/200M_0.0100_AB.dat")
latBB2 = read_lattice("data/0.0195706_0.01/200M_0.0100_BB.dat")

latAA3 = read_lattice("data/0.0189290/200M_0.0050_AA.dat")
latBA3 = read_lattice("data/0.0189290/200M_0.0050_BA.dat")
latAB3 = read_lattice("data/0.0189290/200M_0.0050_AB.dat")
latBB3 = read_lattice("data/0.0189290/200M_0.0050_BB.dat")

latAA4 = read_lattice("data/0.0195706/200M_0.0050_AA.dat")
latBA4 = read_lattice("data/0.0195706/200M_0.0050_BA.dat")
latAB4 = read_lattice("data/0.0195706/200M_0.0050_AB.dat")
latBB4 = read_lattice("data/0.0195706/200M_0.0050_BB.dat")


ax1 = subplot(221,aspect=1)
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

# Create the zoomed-in inset
axins1 = ax1.inset_axes([0.4, 0.55, 0.4, 0.4], xlim=(-200,200), ylim=(-200,200))
try axins1.scatter(latAA1[:,1], latAA1[:,2], s=20, color="c")
catch e
    println("AA lattice is empty!")
end
try axins1.scatter(latBA1[:,1], latBA1[:,2], s=20, color="green")
catch e
    println("BA lattice is empty!")
end
try axins1.scatter(latAB1[:,1], latAB1[:,2], s=20, color="orange")
catch e
    println("AB lattice is empty!")
end
try axins1.scatter(latBB1[:,1], latBB1[:,2], s=20, color="red")
catch e
    println("BB lattice is empty!")
end
axins1.yaxis.set_label_position("right")
axins1.yaxis.tick_right()
axins1.tick_params(labelbottom=false, bottom=false)
ax1.indicate_inset_zoom(axins1, edgecolor="black")



ax2 = subplot(222,aspect=1)
try ax2.scatter(latAA2[:,1], latAA2[:,2], s=10, color="c")
catch e
    println("AA lattice is empty!")
end
try ax2.scatter(latBA2[:,1], latBA2[:,2], s=10, color="green")
catch e
    println("BA lattice is empty!")
end
try ax2.scatter(latAB2[:,1], latAB2[:,2], s=10, color="orange")
catch e
    println("AB lattice is empty!")
end
try ax2.scatter(latBB2[:,1], latBB2[:,2], s=10, color="red")
catch e
    println("BB lattice is empty!")
end
ax2.set_xlim([-3500, 3500])
ax2.set_ylim([-3500, 3500])
#legend(["AA", "BA", "AB", "BB"])
ax2.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)

axins2 = ax2.inset_axes([0.4, 0.55, 0.4, 0.4], xlim=(-200,200), ylim=(-200,200))
try axins2.scatter(latAA2[:,1], latAA2[:,2], s=20, color="c")
catch e
    println("AA lattice is empty!")
end
try axins2.scatter(latBA2[:,1], latBA2[:,2], s=20, color="green")
catch e
    println("BA lattice is empty!")
end
try axins2.scatter(latAB2[:,1], latAB2[:,2], s=20, color="orange")
catch e
    println("AB lattice is empty!")
end
try axins2.scatter(latBB2[:,1], latBB2[:,2], s=20, color="red")
catch e
    println("BB lattice is empty!")
end
axins2.yaxis.set_label_position("right")
axins2.yaxis.tick_right()
axins2.tick_params(labelbottom=false, bottom=false)
ax2.indicate_inset_zoom(axins2, edgecolor="black")



ax3 = subplot(223,aspect=1)
try ax3.scatter(latAA3[:,1], latAA3[:,2], s=10, color="c")
catch e
    println("AA lattice is empty!")
end
try ax3.scatter(latBA3[:,1], latBA3[:,2], s=10, color="green")
catch e
    println("BA lattice is empty!")
end
try ax3.scatter(latAB3[:,1], latAB3[:,2], s=10, color="orange")
catch e
    println("AB lattice is empty!")
end
try ax3.scatter(latBB3[:,1], latBB3[:,2], s=10, color="red")
catch e
    println("BB lattice is empty!")
end
ax3.set_xlim([-3500, 3500])
ax3.set_ylim([-3500, 3500])
#legend(["AA", "BA", "AB", "BB"])
ax3.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)

axins3 = ax3.inset_axes([0.4, 0.55, 0.4, 0.4], xlim=(-200,200), ylim=(-200,200))
try axins3.scatter(latAA3[:,1], latAA3[:,2], s=20, color="c")
catch e
    println("AA lattice is empty!")
end
try axins3.scatter(latBA3[:,1], latBA3[:,2], s=20, color="green")
catch e
    println("BA lattice is empty!")
end
try axins3.scatter(latAB3[:,1], latAB3[:,2], s=20, color="orange")
catch e
    println("AB lattice is empty!")
end
try axins3.scatter(latBB3[:,1], latBB3[:,2], s=20, color="red")
catch e
    println("BB lattice is empty!")
end
axins3.yaxis.set_label_position("right")
axins3.yaxis.tick_right()
axins3.tick_params(labelbottom=false, bottom=false)
ax3.indicate_inset_zoom(axins3, edgecolor="black")





ax4 = subplot(224,aspect=1)
try ax4.scatter(latAA4[:,1], latAA4[:,2], s=10, color="c")
catch e
    println("AA lattice is empty!")
end
try ax4.scatter(latBA4[:,1], latBA4[:,2], s=10, color="green")
catch e
    println("BA lattice is empty!")
end
try ax4.scatter(latAB4[:,1], latAB4[:,2], s=10, color="orange")
catch e
    println("AB lattice is empty!")
end
try ax4.scatter(latBB4[:,1], latBB4[:,2], s=10, color="red")
catch e
    println("BB lattice is empty!")
end
ax4.set_xlim([-3500, 3500])
ax4.set_ylim([-3500, 3500])
#legend(["AA", "BA", "AB", "BB"])
ax4.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)

axins4 = ax4.inset_axes([0.4, 0.55, 0.4, 0.4], xlim=(-200,200), ylim=(-200,200))
try axins4.scatter(latAA4[:,1], latAA4[:,2], s=20, color="c")
catch e
    println("AA lattice is empty!")
end
try axins4.scatter(latBA4[:,1], latBA4[:,2], s=20, color="green")
catch e
    println("BA lattice is empty!")
end
try axins4.scatter(latAB4[:,1], latAB4[:,2], s=20, color="orange")
catch e
    println("AB lattice is empty!")
end
try axins4.scatter(latBB4[:,1], latBB4[:,2], s=20, color="red")
catch e
    println("BB lattice is empty!")
end
axins4.yaxis.set_label_position("right")
axins4.yaxis.tick_right()
axins4.tick_params(labelbottom=false, bottom=false)
ax4.indicate_inset_zoom(axins4, edgecolor="black")



ax1.set_title(L"$\theta=1.0845^{\circ}$, $\delta=0.01 \AA$") 
ax2.set_title(L"$\theta=1.1213^{\circ}$, $\delta=0.01 \AA$")
ax3.set_title(L"$\theta=1.0845^{\circ}$, $\delta=0.005 \AA$")
ax4.set_title(L"$\theta=1.1213^{\circ}$, $\delta=0.005 \AA$")

legend(["AA"], loc="center left", bbox_to_anchor=(1, 0.5))

show()

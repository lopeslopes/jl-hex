include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using Printf

n = 160
a = 2.46

lat_angle = pi/3.0
a1 = [a, 0.0]
a2 = [a*cos(lat_angle), a*sin(lat_angle)]

d = sqrt((a^2)/(2.0*(1.0-cos(2.0*lat_angle))))
d1 = [d*cos(lat_angle/2.0), d*sin(lat_angle/2.0)]

latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)
latA2 = zeros(n ÷ 2, 2)
latB2 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)

ax1 = subplot(131,aspect=1)
ax1.scatter(latA1[:,1], latA1[:,2], s=60, color="blue")
ax1.scatter(latB1[:,1], latB1[:,2], s=60, color="orange")

ax1.quiver(0.0, 0.0, a1[1], a1[2], angles="xy", scale_units="xy", scale=1)
ax1.quiver(0.0, 0.0, a2[1], a2[2], angles="xy", scale_units="xy", scale=1)

ax1.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)
ax1.set_xlim([-4.2, 5.4])
ax1.set_ylim([-4.5, 4.5])

ax2 = subplot(132,aspect=1)
ax2.scatter(latA1[:,1], latA1[:,2], s=60, color="blue")
ax2.scatter(latB1[:,1], latB1[:,2], s=60, color="orange")

ax2.quiver(0.0, 0.0, d1[1], d1[2], angles="xy", scale_units="xy", scale=1)

ax2.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)
ax2.set_xlim([-4.2, 5.4])
ax2.set_ylim([-4.5, 4.5])

ax3 = subplot(133,aspect=1)
ax3.scatter(latA1[:,1], latA1[:,2], s=60, color="blue")
ax3.scatter(latB1[:,1], latB1[:,2], s=60, color="orange")

ax3.quiver(0.0, 0.0, a1[1], a1[2], angles="xy", scale_units="xy", scale=1)
ax3.quiver(0.0, 0.0, a2[1], a2[2], angles="xy", scale_units="xy", scale=1)

ax3.tick_params(left=false, right=false, labelleft=false, labelbottom=false, bottom=false)
ax3.set_xlim([-4.2, 5.4])
ax3.set_ylim([-4.5, 4.5])

show()

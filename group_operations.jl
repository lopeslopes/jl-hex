include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Base.Threads


tol = 5.0e-3
name = @sprintf("%6.4f", tol)
origin = [0.0, 0.0]

latAA = read_lattice("data/0.0191434/"*name*"_AA.dat")
latBA = read_lattice("data/0.0191434/"*name*"_BA.dat")
latAB = read_lattice("data/0.0191434/"*name*"_AB.dat")
latBB = read_lattice("data/0.0191434/"*name*"_BB.dat")

treeAB = KDTree(transpose(latAB))
ind, dist = knn(treeAB, origin, 500)
lenBB = size(latBB)[1]
lenBA = size(latBA)[1]

# FINDING POINTS THAT ALIGN
point_AB = [0.0, 0.0]
point_BB = [0.0, 0.0]
point_BA = [0.0, 0.0]
slope_final = 0.0
for i in 1:500
    point = latAB[ind[i],:]
    slope = point[2]/point[1]
    for j in 1:lenBB
        y_slope = latBB[j,1]*slope
        diff = abs(latBB[j,2]-y_slope)
        if diff<1
            println("    diff: ", diff)
            global point_AB = point
            global point_BB = latBB[j,:]
            global slope_final = point_AB[2]/point_AB[1]
            println("Slope of line: ", slope_final)
            # TEST BA POINTS TOO
            for k in 1:lenBA
                y_BA_slope = latBA[k,1]*slope
                diff_BA = abs(latBA[k,2]-y_BA_slope)
                if diff_BA<1
                    println("    diffBA: ", diff_BA)
                    global point_BA = latBA[k,:]
                end
            end
        end
    end
end

# FIGURING OUT AB POINT COORDINATES IN RELATION TO LATTICE VECTORS A1 AND A2
lat_angle = pi/3.0
a = 2.46
a1 = [a, 0.0]
a2 = [a*cos(lat_angle), a*sin(lat_angle)]

d = sqrt((a^2)/(2.0*(1.0-cos(2.0*lat_angle))))
d1 = [d*cos(lat_angle/2.0), d*sin(lat_angle/2.0)]
origin_a = [0.0, 0.0]
origin_b = origin_a + d1

qtd_a2_AB = round((point_AB[2] - origin_a[2])/a2[2])
AB_aux = origin_a + qtd_a2_AB*a2
qtd_a1_AB = round((point_AB[1] - AB_aux[1])/a1[1])
println("AB point: ", point_AB)
println("    Number of a1: ", qtd_a1_AB)
println("    Number of a2: ", qtd_a2_AB)
AB_test = origin_a + qtd_a1_AB*a1 + qtd_a2_AB*a2
println("    Point obtained with lattice vectors: ", AB_test)
dist_AB = sqrt(point_AB[1]^2 + point_AB[2]^2)
println("    Distance from origin: ", dist_AB)

# FIGURING OUT BB POINT COORDINATES IN RELATION TO LATTICE VECTORS A1 AND A2
qtd_a2_BB = round((point_BB[2] - origin_b[2])/a2[2])
BB_aux = origin_b + qtd_a2_BB*a2
qtd_a1_BB = round((point_BB[1] - BB_aux[1])/a1[1])
println("BB point: ", point_BB)
println("    Number of a1: ", qtd_a1_BB)
println("    Number of a2: ", qtd_a2_BB)
BB_test = origin_b + qtd_a1_BB*a1 + qtd_a2_BB*a2
println("    Point obtained with lattice vectors: ", BB_test)
dist_BB = sqrt(point_BB[1]^2 + point_BB[2]^2)
println("    Distance from origin: ", dist_BB)

println("Ratio distAB/distBB: ", (dist_AB+dist_BB)/dist_BB)

# LINE CROSSING AA, AB AND BB POINTS
x_vals = -12500:12500
y_vals = slope_final*x_vals

# PLOTTING ALL AA, BA, AB AND BB POINTS
ax1 = subplot(111, aspect=1)
ax1.scatter(latAA[:,1], latAA[:,2], color=("green", 0.5), s=8)
ax1.scatter(latBA[:,1], latBA[:,2], color=("red", 0.5), s=8)
ax1.scatter(latAB[:,1], latAB[:,2], color=("orange", 0.5), s=8)
ax1.scatter(latBB[:,1], latBB[:,2], color=("blue", 0.5), s=8)

# SELECTED AB AND BB POINTS
ax1.scatter(0.0, 0.0, color="green", s=50)
ax1.scatter(point_AB[1], point_AB[2], color="magenta", s=50)
ax1.scatter(point_BB[1], point_BB[2], color="purple", s=50)

ax1.scatter(point_BA[1], point_BA[2], color="cyan", s=50)

# PLOTTING LINE
ax1.plot(x_vals, y_vals)

# SETTING PLOT LIMITS AND LEGEND
ax1.set_xlim([-12250, 12250])
ax1.set_ylim([-12250, 12250])
legend(["AA", "BA", "AB", "BB", "selected AA", "selected AB", "selected BB", "selected_BA"])

show()

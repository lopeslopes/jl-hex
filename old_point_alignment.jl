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

latAA = read_lattice("data/0.0191434/1B_"*name*"_AA.dat")
latBA = read_lattice("data/0.0191434/1B_"*name*"_BA.dat")
latAB = read_lattice("data/0.0191434/1B_"*name*"_AB.dat")
latBB = read_lattice("data/0.0191434/1B_"*name*"_BB.dat")

lenAB = size(latAB)[1]
lenBA = size(latBA)[1]
lenBB = size(latBB)[1]
lenAA = size(latAA)[1]
treeBB = KDTree(transpose(latBB))
ind, dist = knn(treeBB, origin, lenBB)

# FINDING BB POINTS THAT ALIGN
point_AB = [0.0, 0.0]
point_BB = [0.0, 0.0]
point_BA = [0.0, 0.0]
slope_final = 0.0
tol_line = 0.1
for i in 1:lenBB
    point = latBB[ind[i],:]
    slope = point[2]/point[1]
    for j in 1:lenAB
        y_slope = latAB[j,1]*slope
        diff_AB = abs(latAB[j,2]-y_slope)
        if diff_AB<tol_line
            # TESTING BA POINTS TOO
            for k in 1:lenBA
                y_BA_slope = latBA[k,1]*slope
                diff_BA = abs(latBA[k,2]-y_BA_slope)
                if diff_BA<tol_line
                    println("    diff  : ", diff_AB)
                    println("    diffBA: ", diff_BA)
                    global point_BA = latBA[k,:]
                    global point_BB = point
                    global point_AB = latAB[j,:]
                    global slope_final = point_BB[2]/point_BB[1]
                    println("Slope of line: ", slope_final)
                end
            end
        end
    end
end

println("------------------------------------------")

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

# FIGURING OUT BA POINT COORDINATES IN RELATION TO LATTICE VECTORS A1 AND A2
qtd_a2_BA = round((point_BA[2] - origin_b[2])/a2[2])
BA_aux = origin_b + qtd_a2_BA*a2
qtd_a1_BA = round((point_BA[1] - BA_aux[1])/a1[1])
println("BA point: ", point_BA)
println("    Number of a1: ", qtd_a1_BA)
println("    Number of a2: ", qtd_a2_BA)
BA_test = origin_b + qtd_a1_BA*a1 + qtd_a2_BA*a2
println("    Point obtained with lattice vectors: ", BA_test)
dist_BA = sqrt(point_BA[1]^2 + point_BA[2]^2)
println("    Distance from origin: ", dist_BA)
# LINE CROSSING AA, AB AND BB POINTS
x_vals = -30000:30000
y_vals = slope_final*x_vals

# PLOTTING ALL AA, BA, AB AND BB POINTS
ax1 = subplot(111, aspect=1)
ax1.scatter(latAA[:,1], latAA[:,2], color=("green", 0.3), s=8)
ax1.scatter(latBA[:,1], latBA[:,2], color=("red", 0.3), s=8)
ax1.scatter(latAB[:,1], latAB[:,2], color=("orange", 0.3), s=8)
ax1.scatter(latBB[:,1], latBB[:,2], color=("blue", 0.3), s=8)

# SELECTED AB AND BB POINTS
ax1.scatter(0.0, 0.0, color="green", s=50)
ax1.scatter(point_AB[1], point_AB[2], color="orange", s=50)
ax1.scatter(point_BB[1], point_BB[2], color="blue", s=50)
ax1.scatter(point_BA[1], point_BA[2], color="red", s=50)

# PLOTTING LINE
ax1.plot(x_vals, y_vals)

# SETTING PLOT LIMITS AND LEGEND
ax1.set_xlim([-40000, 40000])
ax1.set_ylim([-40000, 40000])
legend(["AA", "BA", "AB", "BB", "selected AA", "selected AB", "selected BB", "selected_BA"])

show()

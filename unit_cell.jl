include("hex_utils.jl")
using .HexUtils
using LinearAlgebra
using NearestNeighbors


data_path = "data/0.0186237_200k/"

# superlattice vectors determination by article formula
ang_aux = (60.0/180.0) * π
a = 2.46
v1 = [a, 0.0]
v2 = [a*cos(ang_aux), a*sin(ang_aux)]

p, q, j, steps = read_properties(data_path)
angle_i = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
angle_f = acos((3.0*((q-1)^2) - (p^2))/(3.0*((q-1)^2) + (p^2)))
angle = angle_i + (j/steps)*(angle_f  - angle_i)
println("Angle in radians: ", angle)
println("Angle in degrees: ", (angle * 180) / pi)

if (j==0) || (j==10)
    println("Exact magic angle, calculating supercell vectors...")
    delta = 3/gcd(p, 3)
    gamma = gcd(p+3*q, p-3*q)
    println("Delta: ",delta)
    println("Gamma: ", gamma)
    
    if (delta == 1.0)
        t1 = (1/gamma)*(p+3*q)*v1 + (1/gamma)*(p-3*q)*v2
        t2 = -(1/gamma)*(2*p)*v1 + (1/gamma)*(p+3*q)*v2
    elseif(delta == 3.0)
        t1 = (1/gamma)*(p+q)*v1 - (1/gamma)*(2*q)*v2
        t2 = -(1/gamma)*(p-q)*v1 + (1/gamma)*(p+q)*v2 
    end

    println("|t1| = ", LinearAlgebra.norm(t1))
    println("|t2| = ", LinearAlgebra.norm(t2))
end

moire_period = a/(2*sin(angle/2))
println(" D   = ", moire_period)


latAA = read_lattice(data_path*"latticeAA.dat")
latBA = read_lattice(data_path*"latticeBA.dat")
latAB = read_lattice(data_path*"latticeAB.dat")
latBB = read_lattice(data_path*"latticeBB.dat")


origin = [0.0, 0.0]

# AB range search
try tree = KDTree(transpose(latAB))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("AB points with distance ~D: ", range_ids)
catch e
end

# BA range search
try tree = KDTree(transpose(latBA))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("BA points with distance ~D: ", range_ids)
catch e
end

# AA range search
try tree = KDTree(transpose(latAA))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("AA points with distance ~D: ", range_ids)
catch e
end

# BB range search
try tree = KDTree(transpose(latBB))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("BB points with distance ~D: ", range_ids)
catch e
end

######### doubling the moire period just to check ##################
moire_period = 2*moire_period
# AB range search
try tree = KDTree(transpose(latAB))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("AB points with distance ~2*D: ", range_ids)
catch e
end


# BA range search
try tree = KDTree(transpose(latBA))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("BA points with distance ~2*D: ", range_ids)
catch e
end

# AA range search
try tree = KDTree(transpose(latAA))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("AA points with distance ~2*D: ", range_ids)
catch e
end

# BB range search
try tree = KDTree(transpose(latBB))
    range_ids_outer = inrange(tree, origin, moire_period+0.1)
    range_ids_inner = inrange(tree, origin, moire_period-0.1)
    range_ids = setdiff(range_ids_outer, range_ids_inner)
    println("BB points with distance ~2*D: ", range_ids)
catch e
end

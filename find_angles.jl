include("hex_utils.jl")
using .HexUtils
using PyCall
pygui(:qt5)
using PyPlot
using NearestNeighbors
using Printf
using Base.Threads


# INITIAL DEFINITIONS
n = 10000000
a = 2.46
origin = [0.0, 0.0]

# ALLOCATION OF LATTICES AND FIRST CREATION
println("Creating lattices...")
latA1 = zeros(n ÷ 2, 2)
latB1 = zeros(n ÷ 2, 2)
latA2 = zeros(n ÷ 2, 2)
latB2 = zeros(n ÷ 2, 2)

HexUtils.create_honeycomb_lattice!(latA1, latB1, a, false)
HexUtils.create_honeycomb_lattice!(latA2, latB2, a, false)

distsA1 = zeros(n ÷ 2)
distsB1 = zeros(n ÷ 2)
distsA2 = zeros(n ÷ 2)
distsB2 = zeros(n ÷ 2)
for i in 1:div(n,2)
    distsA1[i] = sqrt(latA1[i,1]^2 + latA1[i,2]^2)
    distsA2[i] = sqrt(latA2[i,1]^2 + latA2[i,2]^2)
    distsB1[i] = sqrt(latB1[i,1]^2 + latB1[i,2]^2)
    distsB2[i] = sqrt(latB2[i,1]^2 + latB2[i,2]^2)
end

# TEST SECTION: TREES
treeA1 = KDTree(transpose(distsA1), reorder=false)
treeB1 = KDTree(transpose(distsB1), reorder=false)
 
tol = 5.0e-4
BA = []
AB = []
println("Tolerance:        ", tol)
for i in 1:div(n,2)
    indBA, distBA = nn(treeB1, [distsA2[i]])
    indAB, distAB = nn(treeA1, [distsB2[i]])
    if distBA < tol
        angA = atan(latA2[i,2]/latA2[i,1])
        angB = atan(latB1[indBA,2]/latB1[indBA,1])
        angle = max(angA, angB) - min(angA, angB)
        if 0.01911 < angle < 0.01917
            push!(BA, angle)
        end
    end
    if distAB < tol
        angA = atan(latA1[indAB,2]/latA1[indAB,1])
        angB = atan(latB2[i,2]/latB2[i,1])
        angle = max(angA, angB) - min(angA, angB)
        if 0.01911 < angle < 0.01917
            push!(AB, angle)
        end
    end
end

println(BA)
println(AB)

# NEW ANGLES OBTAINED FOR TESTING:
# 0.01913261073837491
# 0.01915472274745922
# 0.01914389065354971
# 0.019155348627185043
# 0.019132147632717533
# 0.01916760530700501
# 0.019156276012114448
# 0.019120090428653258
# 0.01914440883336299
# 0.019112562931286936
# 0.01912540230895443
# 0.019151034345490703
# 0.0191182535607608
# 0.019145028761298644
# 0.019145783635169167
# 0.019122613086881235
# 0.01911391097783821
# 0.019163320578591825
# 0.019165942385689805
# 0.019165942385689916
# 0.0191633205785936
# 0.01911391097783499
# 0.01912261308687535
# 0.019145783635161617
# 0.019145028761285543
# 0.019118253560766574
# 0.019151034345471163
# 0.01912540230893478
# 0.019112562931275612
# 0.019144408833343896
# 0.019120090428629055
# 0.019156276012105344
# 0.019167605306979363
# 0.01913214763273996
# 0.019155348627157176
# 0.01914389065352573
# 0.01915472274748997
# 0.019132610738401334
# 0.019132610738401334
# 0.01915472274748997
# 0.01914389065352573
# 0.019160937714272297
# 0.019167605306979585
# 0.01913214763273996
# 0.019120090428629055
# 0.019144408833343896
# 0.019118253560767018
# 0.0191450287612851
# 0.019124340243658833
# 0.019130939585746853
# 0.019153941607066072
# 0.019163320578590437
# 0.019158534293419738
# 0.019155348627157176
# 0.019165942385688917
# 0.01914345108312343
# 0.019159853318520437

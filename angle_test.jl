include("hex_utils.jl")
using .HexUtils

for q in 50:200
    println(q, " ", (180/pi)*magic_angle(1,q))
end

# q1 = 82
# q2 = 81
# angle1 = magic_angle(1,q1)
# angle2 = magic_angle(1,q2)
# steps = 200
# for i in 1:steps
#     angle = angle1 + (angle2-angle1)*i*(1/steps)
#     println(angle)
# end
#
# m = 0
# n = 250
# println(asin((0.5+m)/(3*n+3)))

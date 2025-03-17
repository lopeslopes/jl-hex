# p = 1.0
# q = 3.0

for p in 1:100
    for q in 1:100
        angle = acos((3.0*(q^2) - (p^2))/(3.0*(q^2) + (p^2)))
        if (angle < 0.02) & (angle > 0.0186)
            println("p = ", p, ", q = ", q)
            println("    Angle in radians: ", angle)
            println("    Angle in degrees: ", (angle * 180) / pi)
        end
    end
end

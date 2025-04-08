using PyCall: LinearAlgebra

angle = (60.0/180.0) * π

a = 2.46
a1 = [a, 0.0]
a2 = [a*cos(angle), a*sin(angle)]

a1_3d = [a1[1], a1[2], 0.0]
a2_3d = [a2[1], a2[2], 0.0]
area_lat_cell = LinearAlgebra.norm(LinearAlgebra.cross(a1_3d, a2_3d))

println(area_lat_cell)

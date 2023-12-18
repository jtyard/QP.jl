# Build the 19abcd fields
#
# D = (19-3)(19+1) = 320 
# D0 = 5
# f = 8
# a - adjoin conductor 8 ring class field
# b - adjoin conductor 4 ring class field
# c - adjoin conductor 2 ring class field
# d is the ray class field
#
S = SicData(19)
K = S.K 
OK = S.OK 

rcf = ray_class_field(64*S.OK,S.inf)

rcg, ph = ray_class_group(64*S.OK,S.inf)

art = artin_map(rcf)

PZ = [ph(g) for g in rcg if is_principal_fac_elem(ph(g))[1]]
length(PZ)

a = PZ[2]

art.map1(a)
# how to put a in rcf???
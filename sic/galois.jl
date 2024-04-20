using Oscar, QP

Pa = fiducial("7a")
Pb = fiducial("7b")

S = SicData(7)

F = S.F 

rcf = S.rcf

mu = dwork_modulus(Pb)

art = artin_map(rcf)

art1 = art.map1
art2 = art.map2

gal = domain(art2)

#stab = sub(gal,[g for g in gal if art2(g)(mu) == mu])[1]

#ng, m = norm_group(rcf)

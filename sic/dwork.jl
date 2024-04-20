using Oscar, QP

d = 4
S = SicData(d)

Psi = fiducial(d)

F = S.F 

rcf = S.rcf

mu = dwork_modulus(Psi)

art = artin_map(rcf)

art1 = art.map1
art2 = art.map2

gal = domain(art2)

stab = sub(gal,[g for g in gal if art2(g)(mu) == mu])[1]

ng, m = norm_group(rcf)

hcf = fixed_field(rcf,sub(ng,[g for g in ng if art(m(g))(mu) == mu])[1])
using Oscar, QP


# Checking the modulus of the Dwork hypersurface on which the fiducial lives. 


# Match the Artin map with the Weil representation
to_abs(F) = inv(absolute_simple_field(F)[2])

S = SicData(5)
K = S.K
F = S.F
z5 = zetaN(5,F)

OK = S.OK

P = fiducial(5)

mu = dwork_modulus(P[1,:])

art = artin_map(S.rcf)

G = domain(art.map2)
gal(g) = art.map2(G(g))

println([g for g in G if art.map2(g)(mu) == mu])

f = minpoly(to_abs(F)(mu))
L = number_field(f)[1]
LL = compositum(quadratic_field(5)[1],quadratic_field(-1)[1])[1]
disc = absolute_discriminant(maximal_order(L))

disc, factor(disc), is_isomorphic(L,LL)
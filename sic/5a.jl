using Oscar, QP

# Match the Artin map with the Weil representation


S = SicData(5)
K = S.K
OK = S.OK
s = gen(S.K)

P = fiducial(5)

mu = dwork_modulus(P[1,:])


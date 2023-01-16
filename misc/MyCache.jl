module MyCache

using Oscar

export Z2, ZN

ZN(N) = ResidueRing(ZZ,N)
Z2 = ZN(2)

end
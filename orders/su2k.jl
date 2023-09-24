#################
# Algebras and orders for SU(2)_k
#################
export Dn, su2levelk

# Some orders from https://arxiv.org/abs/1504.04350 (Section 4) and https://arxiv.org/abs/1510.03888 for SU(2)_k

# Dn(n) is totally positive and sqrt(-Dn(n)) generates Q(zeta_n):
#Dn(n) = n % 4 == 0 ? 1 : 4-CyclotomicRealSubfield(n)[2]^2

# q = zeta_{k+2}, q^{1/2} = zeta_{2k + 4}

# julia> is_principal(su2levelk(4-2).D)
# (true, 4)

is_principal(su2levelk(5-2).D)
# (true, 5)

# julia> is_principal(su2levelk(6-2).D)
# (true, 4)

# julia> is_principal(su2levelk(7-2).D)
# (true, 1)

# julia> is_principal(su2levelk(8-2).D)
# (true, 2)


#su2levelkD(k) = discriminant(maximal_order(su2levelk(k)))

#su2levelkR(k) = base_ring(maximal_order(su2levelk(k)))
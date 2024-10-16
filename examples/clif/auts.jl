using Oscar, QP

D8 = dihedral_group(8)
Q8 = quaternion_group(8)

F,i = cyclotomic_field(4)
L = hermitian_lattice(F; gram = identity_matrix(F,2)) 

K,s = quadratic_field(2)
OK = maximal_order(K)
M = quadratic_lattice(K; gram = identity_matrix(K,2)) 


#order(automorphism_group(L))
# 32


#O4 = automorphism_group(integer_lattice(;gram = identity_matrix(QQ,4)))

#order(O4)
# 384


#order(automorphism_group(hermitian_lattice(F; gram = identity_matrix(F,4))))
#6144


#M  = lattice(OK;gram = matrix(OK,2,2,[2,s,s,2])) # doesn't work
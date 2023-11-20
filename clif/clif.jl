using Oscar
using QP 

# Experiments with constructing the 2-qubit Clifford group and its subgroups that stabilize quadratic forms over F_2.

# (Z/2)^4 -> PC_2 -> Sp(4,2) # 720*16 = 11520

# ⟨ζ_8,X_i,Z_i⟩ -> C_2 -> Sp(4,2) # 720*16*8 = 92160

println(order(Sp(4,2))) # 720 elements preserving the bilinear form b = [0 1; 1 0]^(⊕ 2)

# more than one q satisfies b = q + q^T (freedom on the diagonal) qij(x0,x1) = i x0 + j x1 + x0 x1
# O(q00) ≃ O(q01) ≃ O(q10) ≃ Z/3

# What are the right types for finite fields?
# gfp_mat, gfp_elem from GF(2)
# fq_nmod_mat, fq_nmod from GF(2,1)

gf2 = GF(2,1)
q00 = gf2[0 1; 0 0]
q11 = gf2[1 1; 0 1]
o = 0*q00
q = [q00 o; o q00]
b = invariant_bilinear_forms(Sp(4,2))[1]

println(length([g for g in Sp(2,2) if quadratic_form(q11) == quadratic_form(q11)^g]))
println(length([g for g in Sp(2,2) if quadratic_form(q00) == quadratic_form(q00)^g]))

println(length([g for g in Sp(4,2) if transpose(g.elm)*b*g.elm == b]))

#println(length([g for g in Sp(4,2) if quadratic_form( [q00 o; o q00]) == quadratic_form( [q00 o; o q00])^(W*g*transpose(W))]))
#println(length([g for g in Sp(4,2) if quadratic_form( [q00 o; o q11]) == quadratic_form( [q00 o; o q11])^(W*g*transpose(W))]))

# Here is some of the way towards defining affine extension (Z/2)^4 rtimes Sp(4,2) by isomorphism to abelian group

A = abelian_group([2,2,2,2])
V = GF(2,1)^4 #F^n as matrix space of column vectors over F is defined in QP.jl should I submit it to Oscar?

AtoV(a) = V([GF(2,1)(a[i]) for i in 1:4])
VtoA(v) = A([v[i] == GF(2,1)(0) ? 0 : 1 for i in 1:4])

# But maybe a better way is just to define a new matrix group: see affine_group(G) in src/Groups.jl

println(order(affine_group(Sp(4,2))))
println(2^4*720)

    
D = dihedral_group(8)
C,ph = center(D)
c = ph(C[1]) # generator of the center of D8
println(length([s for s in aut(D) if s(c) == c]))

Q = quaternion_group(8)
Z,phh = center(Q)
z = phh(Z[1]) # generator of the center of Q
println(length([s for s in aut(Q) if s(z) == z]))



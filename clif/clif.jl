# (Z/2)^4 -> PC_2 -> Sp(4,2) # 720*16 = 11520

# ⟨ζ_8,X_i,Z_i⟩ -> C_2 -> Sp(4,2) # 720*16*8 = 92160

println(order(Sp(4,2))) # 720 elements preserving the bilinear form b = [0 1; 1 0]^(⊕ 2)

# more than one q satisfies b = q + q^T (freedom on the diagonal) qij(x0,x1) = i x0 + j x1 + x0 x1
# O(q00) ≃ O(q01) ≃ O(q10) ≃

q00 = GF(2)[0 1; 0 0]
q11 = GF(2)[1 1; 0 1]
o = 0*q00
q = [q11 o; o q11]




# q-deformed Schur transform for Uq(sl_2) at roots of unity.

# First need to define maps V(2j-1) -> V(2j) ⊗ V(1)

export ketjm 

function ketjm(j,m)
    J = Int(2*j + 1)
    M = Int(2*m + 1)
    [i == M ? 1 : 0 for i in 1:J]
end


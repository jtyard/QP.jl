# q-deformed Schur transform for Uq(sl_2) at roots of unity.

# First need to define maps V(2j-1) -> V(2j) ⊗ V(1)

export ketjm, UqMod, j_plus_halfs, j_minus_halfs, qint

function ketjm(j,m)
    [k == m ? 1 : 0 for k in j:-1:-j]
end

# Maybe rewrite this to live in the real cyclotomic field
# Returns [m] at q^{1/2} = 2nth root of unity
function qint(n,m) 
    if m==0 
        return 0
    end
    q = iseven(n) ? zetaN(2*n)  : -zetaN(n)
    sum([q^(m-1 - 2*i) for i in 0:m-1])
end






struct UqMod
    j
end


# ∣j1+1/2,m⟩∣1/2,±1/2⟩

function j_plus_halfs(j1,m,k) 
    vec = zeros(4 * j1 + 2)
    q4 = zetaN(4*k+4)

    if m < j1
        vec[1 + 2*Int(j1 - (m + 1/2)) + 1] = sqrt(q^(-(j1 + 1/2 + m)/2) * qint(j1 + 1/2 - m,k+2)//qint(2*j1+1,k+2))
    end
    if m > -j1
        vec[1 + 2*Int(j1 - (m - 1/2))] = sqrt(q^((j1 + 1/2 - m)/2) * qint(j1 + 1/2 + m,k+2)//qint(2*j1+1,k+2))
    end 
    vec
end

function j_minus_halfs(j1,m,k)
    vec = zeros(4 * j1 + 2)
    if m < j1
        vec[1 + 2*Int(j1 - (m + 1/2)) + 1] = -sqrt(q^((j1 + 1/2 - m)/2) * qint(j1 + 1/2 + m,k+2)//qint(2*j1+1,k+2))
    end
    if m > -j1
        vec[1 + 2*Int(j1 - (m - 1/2))] = sqrt(q^(-(j1 + 1/2 + m)/2) * qint(j1 + 1/2 + m,k+2)//qint(2*j1+1,k+2))
    end 
    vec
end




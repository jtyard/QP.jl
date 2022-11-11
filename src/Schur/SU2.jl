# q-deformed Schur transform for Uq(sl_2) at roots of unity.

# First need to define maps V(2j-1) -> V(2j) ⊗ V(1)

import Base.show
export Spin, ket,ketjm, UqMod, j_plus_halfs, j_minus_halfs, qint, show

struct Spin
    n::fmpz #units of 1/2 (really hbar omega/ 2) 
    j::fmpq
    Spin(j) = new(fmpz(2*j),fmpz(2*j)//2)
end

function Base.show(io::IO, s::Spin)
    print(io, "Spin-", s.j)
end

struct SpinKet
    j::Spin
    m::Spin
end

function Base.show(io::IO, s::SpinKet)
    print(io, "|", s.j.j,",",s.m.j,"⟩")
end

function ket(j::Spin,m::Spin) 
    SpinKet(j,m)
    #fmpz[Spin(k) == m ? 1 : 0 for k in j.j:-1:-j.j]
end

# The above is a nice idea but it does't really work as expected because e.g. Spin(1) != Spin(1) (== on a struct uses === on the fields for some reason).

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

# Another approach (used to be in Fmatrix.jl)
#_, q = RationalFunctionField(QQ,"q")
#qint(n) = sum([q^(n-1 - 2*i) for i in 0:n-1])





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




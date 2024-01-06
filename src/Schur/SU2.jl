# Experimental - not much useful here yet

# Eventually - q-deformed Schur transform for Uq(sl_2) at roots of unity.

# First need to define maps V(2j-1) -> V(2j) ⊗ V(1)

using Memoize

import Base.show
export Spin, ket,ketjm, UqMod, j_plus_halfs, j_minus_halfs, show

struct Spin
    n::ZZRingElem #units of 1/2 (really ħω/2) 
    j::QQFieldElem
    @memoize Spin(j) = new(ZZRingElem(2*j),ZZRingElem(2*j)//2) # so Spin(1) == Spin(1)
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
    #ZZRingElem[Spin(k) == m ? 1 : 0 for k in j.j:-1:-j.j]
end


function ketjm(j,m)
    [k == m ? 1 : 0 for k in j:-1:-j]
end


struct UqMod
    j
end


# ∣j1+1/2,m⟩∣1/2,±1/2⟩

function j_plus_halfs(j1,m,k) 
    vec = zeros(4 * j1 + 2)
    qq = zetaN(Int(4*k+4))

    if m < j1
        vec[1 + 2*Int(j1 - (m + 1/2)) + 1] = sqrt(qq^(-(2*j1 + 1 + 2*m)) * qint(j1 + 1/2 - m,k+2)//qint(2*j1+1,k+2))
    end
    if m > -j1
        vec[1 + 2*Int(j1 - (m - 1/2))] = sqrt(qq^((2*j1 + 1 - 2*m)) * qint(j1 + 1/2 + m,k+2)//qint(2*j1+1,k+2))
    end 
    vec
end

function j_minus_halfs(j1,m,k)
    vec = zeros(4 * j1 + 2)
    qq = zetaN(Int(4*k+4))

    if m < j1
        vec[1 + 2*Int(j1 - (m + 1/2)) + 1] = -sqrt(qq^((2*j1 + 1 - 2*m)) * qint(Int(j1 + 1/2 + m),k+2)//qint(2*j1+1,k+2))
    end
    if m > -j1
        vec[1 + 2*Int(j1 - (m - 1/2))] = sqrt(qq^(-(2*j1 + 1 + 2*m)) * qint(Int(j1 + 1/2 + m),k+2)//qint(2*j1+1,k+2))
    end 
    vec
end




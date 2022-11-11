# Useful rings, fields and elements
using Oscar

export ZN, zetaN, dagger, complex_conjugate, primes, myorder, Z2, Z3, Z4, Z5, Z6, Z7, Z8

#################
# Ring   Element type
# ZZ     fmpz 
# QQ     fmpq
# ZN(N)  nmod 
#
# Note that characteristic(base_ring(___)) works for GF(p)[ ] and ResidueClassRing(ZZ,4)
#################

Base.Int(a::nmod) = Int(ZZ(a))
ZN(N) = ResidueRing(ZZ,N) 
Z2 = ZN(2)
Z3 = ZN(3)
Z4 = ZN(4)
Z5 = ZN(5)
Z6 = ZN(6)
Z7 = ZN(7)
Z8 = ZN(8)
zetaN(N) = cyclotomic_field(N)[2]

roots_of_unity(F,n) = roots(PolynomialRing(F)[2]^n-1)
primitive_roots_of_unity(F,n) = roots(PolynomialRing(F)[1](collect(coefficients(cyclotomic_field(n)[1].pol))))

zetaN(n,F) = primitive_roots_of_unity(F,n)[1]

# Create row vectors over R (or do I want column vectors??)
import Base.^
R^n = MatrixSpace(R,1,n)

# note that rand(R^n) actually works

# conjugate_transpose

dagger(M) = transpose( map(complex_conjugation(base_ring(M)), M ) )

complex_conjugate(x) = complex_conjugation(parent(x))(x)


import Oscar.primes

primes(N::Int) = [n for n in 2:N if is_prime(n)]

primes(M::Int,N::Int) = [n for n in M:N if is_prime(n)]/9

function myorder(a)
    aa = a
    id = a^0
    for i in 1:1000
        if aa == id
            return i
        end
        aa = aa*a
    end
    return Inf
end 


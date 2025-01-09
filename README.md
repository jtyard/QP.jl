# QP.jl - Quantum Programming

Experimental integrated open-source Julia/Oscar/Nemo/Hecke/GAP/Polymake/ANTIC workflow for quantum computing research as alternative to Python/Magma/Sage/Pari.  Heavily using [`Oscar.jl`](https://docs.oscar-system.org/stable/) which in turn wraps [`GAP.jl`](https://github.com/oscar-system/GAP.jl) for group theory, [`Polymake.jl`](https://github.com/oscar-system/Polymake.jl) for polyhedral geometry and [`Singular.jl`](https://github.com/oscar-system/Singular.jl) for algebraic geometry and invariant theory.  Oscar incorporates [`Hecke.jl`](https://github.com/thofma/Hecke.jl) for computational algebraic number theory and class field theory, which wraps [`ANTIC`](https://github.com/flintlib/antic) for fast number theory computations in C.



```julia
julia> using Oscar, QP

julia> ZN(3)
Integers modulo 3

julia> Z3
Integers modulo 3

julia> Z3[4 2]
[1   2]

julia> ket(Z3[0 0]) + ket(Z3[1 1]) + ket(Z3[2 2])
[1]
[0]
[0]
[0]
[1]
[0]
[0]
[0]
[1]

julia> heis(Z3[4 2])
[  0   0   -z_3 - 1]
[z_3   0          0]
[  0   1          0]

julia> weil_U(Z5[2 1; 0 3])
[-1                         0      0      0                         0]
[ 0                         0   -z_5      0                         0]
[ 0                         0      0      0   z_5^3 + z_5^2 + z_5 + 1]
[ 0   z_5^3 + z_5^2 + z_5 + 1      0      0                         0]
[ 0                         0      0   -z_5                         0]

julia> X = gen(matrix_polynomial_ring(QQ,8))
[X_{0,0}   X_{0,1}   X_{0,2}   X_{0,3}   X_{0,4}   X_{0,5}   X_{0,6}   X_{0,7}]
[X_{1,0}   X_{1,1}   X_{1,2}   X_{1,3}   X_{1,4}   X_{1,5}   X_{1,6}   X_{1,7}]
[X_{2,0}   X_{2,1}   X_{2,2}   X_{2,3}   X_{2,4}   X_{2,5}   X_{2,6}   X_{2,7}]
[X_{3,0}   X_{3,1}   X_{3,2}   X_{3,3}   X_{3,4}   X_{3,5}   X_{3,6}   X_{3,7}]
[X_{4,0}   X_{4,1}   X_{4,2}   X_{4,3}   X_{4,4}   X_{4,5}   X_{4,6}   X_{4,7}]
[X_{5,0}   X_{5,1}   X_{5,2}   X_{5,3}   X_{5,4}   X_{5,5}   X_{5,6}   X_{5,7}]
[X_{6,0}   X_{6,1}   X_{6,2}   X_{6,3}   X_{6,4}   X_{6,5}   X_{6,6}   X_{6,7}]
[X_{7,0}   X_{7,1}   X_{7,2}   X_{7,3}   X_{7,4}   X_{7,5}   X_{7,6}   X_{7,7}]

julia> tr(X^2)
X_{0,0}^2 + 2*X_{0,1}*X_{1,0} + 2*X_{0,2}*X_{2,0} + 2*X_{0,3}*X_{3,0} + 2*X_{0,4}*X_{4,0} + 2*X_{0,5}*X_{5,0} + 2*X_{0,6}*X_{6,0} + 2*X_{0,7}*X_{7,0} + X_{1,1}^2 + 2*X_{1,2}*X_{2,1} + 2*X_{1,3}*X_{3,1} + 2*X_{1,4}*X_{4,1} + 2*X_{1,5}*X_{5,1} + 2*X_{1,6}*X_{6,1} + 2*X_{1,7}*X_{7,1} + X_{2,2}^2 + 2*X_{2,3}*X_{3,2} + 2*X_{2,4}*X_{4,2} + 2*X_{2,5}*X_{5,2} + 2*X_{2,6}*X_{6,2} + 2*X_{2,7}*X_{7,2} + X_{3,3}^2 + 2*X_{3,4}*X_{4,3} + 2*X_{3,5}*X_{5,3} + 2*X_{3,6}*X_{6,3} + 2*X_{3,7}*X_{7,3} + X_{4,4}^2 + 2*X_{4,5}*X_{5,4} + 2*X_{4,6}*X_{6,4} + 2*X_{4,7}*X_{7,4} + X_{5,5}^2 + 2*X_{5,6}*X_{6,5} + 2*X_{5,7}*X_{7,5} + X_{6,6}^2 + 2*X_{6,7}*X_{7,6} + X_{7,7}^2

julia> partial_trace(X,[2,2,2],[1,3])
[X_{0,0} + X_{1,1} + X_{4,4} + X_{5,5}   X_{0,2} + X_{1,3} + X_{4,6} + X_{5,7}]
[X_{2,0} + X_{3,1} + X_{6,4} + X_{7,5}   X_{2,2} + X_{3,3} + X_{6,6} + X_{7,7}]

julia> reduced_operator(X,[2,2,2],[2])
[X_{0,0} + X_{1,1} + X_{4,4} + X_{5,5}   X_{0,2} + X_{1,3} + X_{4,6} + X_{5,7}]
[X_{2,0} + X_{3,1} + X_{6,4} + X_{7,5}   X_{2,2} + X_{3,3} + X_{6,6} + X_{7,7}]

julia> sic(3)
9-element Vector{AbstractAlgebra.Generic.MatSpaceElem{nf_elem}}:
 [0 0 0; 0 1 -1; 0 -1 1]
 [0 0 0; 0 1 -z_3; 0 z_3+1 1]
 [0 0 0; 0 1 z_3+1; 0 -z_3 1]
 [1 -1 0; -1 1 0; 0 0 0]
 [1 -z_3 0; z_3+1 1 0; 0 0 0]
 [1 z_3+1 0; -z_3 1 0; 0 0 0]
 [1 0 -1; 0 0 0; -1 0 1]
 [1 0 z_3+1; 0 0 0; -z_3 0 1]
 [1 0 -z_3; 0 0 0; z_3+1 0 1]

julia> 2*fiducial("7b")[1,:]
7-element Vector{Hecke.RelNonSimpleNumFieldElem{AbsSimpleNumFieldElem}}:
 2
 -_$1 - sqrt(2) - 1
 -_$1 - sqrt(2) - 1
 _$1 - sqrt(2) - 1
 -_$1 - sqrt(2) - 1
 _$1 - sqrt(2) - 1
 _$1 - sqrt(2) - 1

julia> SicData(5)
Constructing number field
Finding complex conjugation
SicData(5, 12, 12, 1, Polynomial ring in 5×5 variables X_{⋅,⋅} over Rational field, Real quadratic field defined by x^2 - 3, InfPlc[Infinite place corresponding to (Complex embedding corresponding to -1.73 of real quadratic field), Infinite place corresponding to (Complex embedding corresponding to 1.73 of real quadratic field)], Maximal order of Real quadratic field defined by x^2 - 3 
with basis AbsSimpleNumFieldElem[1, sqrt(3)], sqrt(3) + 2, Order of Real quadratic field defined by x^2 - 3
with Z-basis AbsSimpleNumFieldOrderElem[1, -sqrt(3) + 6], Order of Real quadratic field defined by x^2 - 3
with Z-basis AbsSimpleNumFieldOrderElem[1, -sqrt(3) + 2], -sqrt(3) + 2, Class field defined mod (<5, 5>, InfPlc{AbsSimpleNumField, AbsSimpleNumFieldEmbedding}[Infinite place corresponding to (Complex embedding corresponding to -1.73 of real quadratic field), Infinite place corresponding to (Complex embedding corresponding to 1.73 of real quadratic field)]) of structure Z/2 x Z/8, Non-simple number field of degree 16 over real quadratic field, Map: non-simple number field -> non-simple number field, #undef, #undef, #undef)

```

Some julia tips can be found [here](julia).

Supported in part by the NSERC Discovery under Grant No. RGPIN-2018-04742, the NSERC/EU project FoQaCiA under Grant No. ALLRP-569582-21, the [Institute for Quantum Computing](https://uwaterloo.ca/institute-for-quantum-computing/) and the [Perimeter Institute for Theoretical Physics](https://perimeterinstitute.ca/).

## Current 
- Rings of polynomial functions on matrices over number fields 
- Explicit computing in bases (vector, operator, etc)
- Useful matrices and tensor products 
- Multiqudit bases parameterized by elements of finite cyclic rings (~~`nmod`~~`zzModRingElem`), of modules over them (~~`nmod_mat`~~`zzModMatrix`) or of finite abelian groups (`FinGenAbGroupElem`)
- Single-qudit generalized Paulis in all dimensions
- Qudit generalized Clifford = Weil representation in prime dimensions
- Projective linear groups 
- Computing properties of SIC-POVMs 

## SIC-POVMs
- `SicData(d)` (or `SicData(d,build_nf=false)` to skip building the number field)
- `fiducial(d)` and `sic(d)` for `d=2,3,3,4,5,7` as well as `fiducial("7a")`, `fiducial("7b")` etc.
- [Class field thery in Hecke](ClassFields.md)


## TODO 
- [Group actions in GAP/Oscar](Groups.md) 
- [Arithmetic of quantum circuits](Arithmetic.md)
- Characteristic and Wigner functions
- Generalize Weil to composite $N$
  - Sort out even case
  - ~~Fix this [Oscar issue](https://github.com/oscar-system/Oscar.jl/issues/649) to implement~~
  ~~$\mathrm{SL}(2,\mathbb{Z}/N)$ for all $N$.~~





## More to do
- Stabilizer/graph states 
- Parametric models / exponential families
- Solving pentagon and hexagon equations 
- Ross-Selinger and extensions
- Schur transform

## Long-term vision
- Abstract framework for modeling quantum systems
- Implementation-agnostic (native julia, arb, FLINT, ..., TensorFlow)
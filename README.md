# QP.jl - Quantum Programming

Experimental integrated open-source Julia/Oscar/Nemo/Hecke/GAP/Polymake/ANTIC workflow for quantum computing research as alternative to Python/Magma/Sage/Pari.  Heavily using [Oscar.jl](https://docs.oscar-system.org/stable/) which in turn wraps [`GAP.jl`](https://github.com/oscar-system/GAP.jl) for group theory, [`Polymake.jl`](https://github.com/oscar-system/Polymake.jl) for polyhedral geometry and [`Singular.jl`](https://github.com/oscar-system/Singular.jl) for algebraic geometry and invariant theory.  Oscar incorporates [`Hecke.jl`](https://github.com/thofma/Hecke.jl) for computational algebraic number theory and class field theory, which wraps [`ANTIC`](https://github.com/flintlib/antic) for fast number theory computations in C.



```julia
julia> using Oscar, QP

julia> ZN(3)
Integers modulo 3

julia> Z3
Integers modulo 3

julia> Z3[4 2]
[1   2]

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

julia> Xij(3)
[X_{0,0}   X_{0,1}   X_{0,2}]
[X_{1,0}   X_{1,1}   X_{1,2}]
[X_{2,0}   X_{2,1}   X_{2,2}]

julia> Xij(3) + Xij(3)

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

julia> SicData(5)
SicData(5, 12, 12, 1, Real quadratic field defined by x^2 - 3, InfPlc[Infinite place corresponding to (Complex embedding corresponding to -1.73 of real quadratic field defined by x^2 - 3), Infinite place corresponding to (Complex embedding corresponding to 1.73 of real quadratic field defined by x^2 - 3)], Maximal order of Real quadratic field defined by x^2 - 3 
with basis nf_elem[1, sqrt(3)], sqrt(3) + 2, Order of Real quadratic field defined by x^2 - 3
with Z-basis NfOrdElem[1, -sqrt(3) + 6], Order of Real quadratic field defined by x^2 - 3
with Z-basis NfOrdElem[1, -sqrt(3) + 11//2], -sqrt(3) + 11//2, Class field defined mod (<5, 5>, InfPlc{AnticNumberField, NumFieldEmbNfAbs}[Infinite place corresponding to (Complex embedding corresponding to -1.73 of real quadratic field defined by x^2 - 3), Infinite place corresponding to (Complex embedding corresponding to 1.73 of real quadratic field defined by x^2 - 3)]) of structure Abelian group with structure: Z/2 x Z/8)

julia> fiducial(5)

```

Some julia tips can be found [here](julia).

Supported in part by the NSERC Discovery under Grant No. RGPIN-2018-04742, the NSERC/EU project FoQaCiA under Grant No. ALLRP-569582-21, the [Institute for Quantum Computing](https://uwaterloo.ca/institute-for-quantum-computing/) and the [Perimeter Institute for Theoretical Physics](https://perimeterinstitute.ca/).

## Current 
- Rings of polynomial functions on matrices over number fields 
- Explicit computing in bases (vector, operator, etc)
- Useful matrices and tensor products 
- Multiqudit bases parameterized by elements of finite cyclic rings (`nmod`), of modules over them (`nmod_mat`) or of finite abelian groups (`GrpAbFinGenElem`)
- Single-qudit generalized Paulis in all dimensions
- Qudit generalized Clifford = Weil representation in prime dimensions
- Projective linear groups 
- Computing properties of SIC-POVMs 

## SIC-POVMs
- `SicData(d)` (or `SicData(d,build_nf=false)` to skip building the number field)
- `fiducial(d)` and `sic(d)` for `d=2,3,3,4,5`


## Class fields
http://www.thofma.com/Hecke.jl/v0.6.1/class_fields/intro.html

Some relevant functions from Hecke:
- Relative automorphism generators from  `Hecke.automorphism_groupQQ`
- Defining ray class fields: `rcf = ray_class_field((iseven(d) ? 2*d : d)*OK,infinite_places(K))` for `K = quadratic_field((d-3)*(d+1))[1]` and `OK = maximal_order(K)`.
- `number_field(rcf)` computes generators `absolute_automorphism_group(rcf)` of the full automorphism group and now works for all `d`.   
- TO DO: check that we can construct the ring ray class fields.
- `MapClassGrp` : {quotient of the class group} -> {ideals} 
- `MapRayClassGrp` : {quotient of a ray class group} -> {ideals prime to the conductor}
- `ClassField` 
- `ClassField_pp` Cyclic class field of prime-power degree
- `artin_map(rcf)` gives map from the ideal group to the set of automorphisms of `number_field(rcf)` i.e. from a `FacElemMon{Hecke.NfAbsOrdIdlSet{AnticNumberField, nf_elem}}` to a `Hecke.NfMorSet{NfRelNS{nf_elem}}`
- `complex_conjugation(F,infplace)` extends the Artin map to infinite places, giving complex conjugation in the corresponding complex embeddings. 
- `automorphism_group(rcf)` gives a map from a `GrpGen` to the set of automorphisms of `number_field(rcf)` fixing the base, and `inv(rcf.quotientmap)` works. 
## TODO 
- Group actions
- Characteristic and Wigner functions
- Generalize Weil to composite $N$
  - Sort out even case
  - Fix this [Oscar issue](https://github.com/oscar-system/Oscar.jl/issues/649) to implement 
  $\mathrm{SL}(2,\mathbb{Z}/N)$ for all $N$.
- Projective schemes - much more development on this by now.

### Towards arithmetic of quantum circuits
Let $R$ be an integral domain and $B$ a central simple algebra over $\mathrm{Frac}(R)$.  
- An **$R$-lattice** is a finitely generated torsion-free $R$-module.
- A **fractional $R$-ideal** $I$ of $B$ is a full-rank $R$-sublattice.
  - **Left order** $`O_\ell(I) = \{ b \in B : bI \subset I \}`$
  - **Right order**  $`O_r(I) = \{ b \in B : bI \subset I \}`$
  - $I$ is an $O_\ell(I)-O_r(I)$-bimodule, a left fractional $O_\ell(I)$ 
ideal and a right fractional $O_r(I)$-ideal.
  - There is a **category** whose objects are the orders and whose morphisms are the fractional ideals, with composition given by proper products.  Fractional ideals **connect** their left/right orders.
- $O\simeq O'$ iff $bOb^{-1} = O'$ for some $b \in B^\times$ 
iff they are connected by a principal ideal $bO = O'b$. So the isomorphism classes are the conjugacy classes.
- An **invertible ideal** of $B$ is a fractional ideal $I$ with an **inverse** $I^{-1}$, satisfying 
$I^{-1}I = O_R(I)$ and therefore $I I^{-1} = O_\ell(I)$.
  - The orders and the invertible ideals form a groupoid called the **core** of the above category.
- A **normal ideal** of $B$ is a fractional $R$-ideal with maximal left/right orders.
  - normal $\Rightarrow$ invertible.  
  -  The maximal orders and normal ideals of $B$ form the **Brandt groupoid**.
- An **integral ideal** of $B$ is a normal $R$-ideal $I$ with $I \subset O_\ell(I) \cap O_r(I)$. 
- An **Eichler order** is the intersection of two maximal orders. By the above definitions, integral ideals are Eichler.  
- A **maximal left ideal** $I$ of $B$ is a maximal left ideal of $O_\ell(I)$.  
- Each integral ideal $I$ of $B$ can be factored into a **proper product** $I = M_1 M_2 \cdots M_n$ of maximal left ideals $M_i$ of $B$ with $O_r(M_i) = O_\ell(M_{i+1})$ for $1 \leq i < n$, 
where $n$ is the length of $O_\ell(I)/I$.




We need the following basic operations to compute with maximal or Eichler quaternion orders (see https://arxiv.org/abs/0808.3833):

- Check isomorphism of fractional ideals: `is_isomorphic` (Oscar) `IsIsomorphic` (Magma), reduces to `is_principal` / `IsPrincipal`
  - Indefinite case: check image in ray class group mod the infinite ramified primes of $B$. 
  - Definite case: Solve shortest lattice vector problem  
- Compute connecting fractional ideals `I(OO,O)` such that left ideal is OO and right ideal is O.

The main difficult tasks are the following:
- Compute representatives for the conjugacy classes = types = isomorphism classes of orders
  - `ConjugacyClasses(O)` in Magma
- Compute representatives for the 2-sided ideal class group 
  - Extends the class group of the base maximal order by square roots of ramified primes
  - `TwoSidedIdealClassGroup(O)` in Magma
- Combining the previous two (KV 2.10) computes representatives 
  `[J*I(OO,O) : OO in ConjugacyClasses(O), J in TwoSidedIdealClassGroup(OO)]` 
  for the set $\mathrm{Pic}_\ell(O)$ of left-equivalence classes of invertible right $O$-ideals. 
  - `RightIdealClasses(O)` in Magma
  - Cardinality is the class number

- Hecke does one split prime.  Can we use `PolyMake.jl` for totally positive and more general? 


## More to do
- Stabilizer/graph states 
- Parametric models / exponential families
- Solving pentagon and hexagon equations 
- Ross-Selinger and extensions
- Schur transform

## Long-term vision
- Abstract framework for modeling quantum systems
- Implementation-agnostic (native julia, arb, FLINT, ..., TensorFlow)
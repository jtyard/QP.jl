# QP.jl - Quantum Programming

Experimental integrated open-source Julia/Oscar/Nemo/Hecke/GAP/Polymake/ANTIC workflow for quantum computing research as alternative to Python/Magma/Sage/Pari.  Heavily using [Oscar.jl](https://docs.oscar-system.org/stable/) which in turn wraps [`GAP.jl`](https://github.com/oscar-system/GAP.jl) for group theory, [`Polymake.jl`](https://github.com/oscar-system/Polymake.jl) for polyhedral geometry and [`Singular.jl`](https://github.com/oscar-system/Singular.jl) for algebraic geometry and invariant theory.  It also incorporates [`Hecke.jl`](https://github.com/thofma/Hecke.jl) for computational algebraic number theory...



```julia
julia> using QP
julia> using Oscar


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
```

Some julia tips can be found [here](julia)

## Current 
- Rings of polynomial functions on matrices over number fields 
- Explicit computing in bases (vector, operator, etc)
- Useful matrices and tensor products 
- Qudit generalized Paulis
- Prime qudit generalized Clifford = Weil representation
- Projective linear groups 
- Computing properties of SIC-POVMs 

## SIC-POVMs
- `SicData(d)` (or `SicData(d,build_nf=true)` to build the number field)



## Class fields
- http://www.thofma.com/Hecke.jl/v0.6.1/class_fields/intro.html
- Relative automorphism generators from  `Hecke.automorphism_groupQQ`
- Things used to go wrong in computing the full automorphism group in many cases e.g. `ray_class_field(5*OK,infinite_places(K))` for `K = quadratic_field(3)` runs forever but they were nice enough to fix for all the ray class fields of interest.  Though I haven't checked the ring class fields and probably won't any time soon.
- `MapClassGrp` : {quotient of the class group} -> {ideals} 
- `MapRayClassGrp` : {quotient of a ray class group} -> {ideals prime to the conductor}
- `ClassField` 
- `ClassField_pp` Cyclic class field of prime-power degree

Given ray class field like `rcf = ray_class_field(5*OK,inf)` and the corresponding relative number field `A = number_field(rcf)` (henceforth is accessible via `rcf.A`), how do I compute the Artin map on ideals and real places? `automorphism_group(rcf)` gives me a map from a `GrpGen` to the set of automorphisms of `A` fixing the base, and `inv(rcf.quotientmap)` (this was fixed and now `artin_map` and `complex_conjugation(F,infplace)` work)
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
- A **fractional $R$-ideal*** $I$ of $B$ is a full-rank $R$-sublattice.
  - **Left order** $`O_L(I) = \{ b \in B : bI \subset I \}`$
  - **Right order**  $`O_L(I) = \{ b \in B : bI \subset I \}`$
  - $I$ is an $O_L(I)-O_R(I)$-bimodule, a left fractional $O_L(I)$ 
ideal and a right fractional $O_R(I)$-ideal.
  - A **proper product** $I_1 I_2 \dotsc I_n$ of fractional ideals has $O_R(I_i) = O_L(I_{i+1})$ for $1 \leq i < n$.
  - There is a category whose objects are the orders and whose morphisms are the fractional ideals, with composition given by proper products.  Fractional ideals **connect** their left/right orders.
- $O\simeq O'$ iff $bOb^{-1} = O'$ for some $b \in B^\times$ 
iff they are connected by a principal ideal $bO = O'b$. So the isomorphism classes are the conjugacy classes.
- An **invertible ideal** is a fractional ideals $I$ with an *inverse* $I^{-1}$, satisfying 
$I^{-1}I = O_R(I)$ therefore $I I^{-1} = O_L(I)$.
  - The orders and the invertible ideals form a groupoid.

- A **normal $R$-ideal** is a fractional $R$-ideal with maximal left/right orders.
  - Every normal ideal is invertible.  
  -  The maximal orders and normal ideals form the **Brandt groupoid**.
- An **integral $R$-ideal** is a normal ideal $I$ with $I \subset O_L(I) \cap O_R(I)$ though it might make sense to drop the normal requirement.  




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
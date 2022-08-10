# QP.jl - Quantum Programming

```julia
julia> using QP
julia> using Oscar


julia> ZN(3)
Integers modulo 3

julia> ZN(3)[4 2]
[1   2]

julia> heis(ZN(3)[4 2])
[  0   0   -z_3 - 1]
[z_3   0          0]
[  0   1          0]

julia> Xij(3)
[X_{0,0}   X_{0,1}   X_{0,2}]
[X_{1,0}   X_{1,1}   X_{1,2}]
[X_{2,0}   X_{2,1}   X_{2,2}]
```

- Abstract framework for modeling quantum systems
- Explicit computing in bases (vector, operator, etc)
- Implementation-agnostic (native julia, arb, FLINT, ..., TensorFlow)
- Julia/Oscar/Nemo/Hecke/ANTIC/GAP/Polymake vs Python/Magma/Sage/Pari

## Current 

- Rings of polynomial functions on matrices over number fields 
- Useful matrices and tensor products 
- Qudit generalized Paulis
- Prime qudit generalized Clifford = Weil representation
- Projective linear groups 
- Investigating properties of SIC-POVMs 

## TODO 
- Group actions
- Projective schemes
- Generalize Weil to composite $N$
  - Sort out even case
  - Fix this [Oscar issue](https://github.com/oscar-system/Oscar.jl/issues/649) to implement 
  $\mathrm{SL}(2,\mathbb{Z}/N)$ for all $N$.



## Towards arithmetic of quantum circuits
Let $R$ be an integral domain and $B$ 
a central simple algebra over Frac$(R)$.  
- An $R$-lattice is a finitely generated torsion-free $R$-module.
- A fractional $R$-ideal of $B$ is a full-rank $R$-sublattice.
Each fractional $R$-ideal $I$ is an $O_L(I)-O_R(I)$-bimodule.  Proper products of fractional ideals have matching left/right orders).  Connected orders are isomorphic = conjugate iff they are connected by a principal ideal. So the conjugacy classes of orders are the isomorphism classes.

- An invertible $R$-ideal $I$ of $B$ is a fractional $R$-ideal $I$ such that there exists a fractional $R$-ideal $J$ with $IJ$  

The invertible fractional ideals form a groupoid.

- An normal $R$-ideal is a fractional $R$-ideal whose left/right orders are maximal.

The normal ideals form the Brandt groupoid.

- An integral $R$-ideal is a normal ideal contained in its left/right orders.

 
The invertible $R$-ideals in a central simple algebra over Frac$(R)$ form a groupoid under proper products (meaning the left and right orders coincide).  The Brandt groupoid consists of the maximal orders 

Normal ideals have O_L and O_R maximal.

Need to be able to do the following basic operations (see Kirchmer & Voight - Algorithmic enumeration of ideal classes for quaternion orders):

- Check isomorphism of fractional ideals: `is_isomorphic` (Oscar) `IsIsomorphic` (Magma) reduces to `is_principal` / `IsPrincipal`
- Compute connecting fractional ideals `I(O,OO)` such that left ideal is O and right ideal is OO.

The main difficult tasks are the following:
- Compute representatives for the conjugacy classes = types = isomorphism classes of orders
  - `ConjugacyClasses(OO)` in Magma
- Compute representatives for the 2-sided ideal class group 
  - `TwoSidedIdealClassGroup(OO)` in Magma
  - Extends the class group of the base maximal order by square roots of ramified primes
- Combining the previous two (KV 2.10) lets us compute representatives `[J*I(O,OO) : O in ConjugacyClasses(OO), J in TwoSidedIdealClassGroup(O)]` for the set $\mathrm{Pic}_\ell(\mathcal{O})$ of left-equivalence classes of invertible right $\mathcal{O}$-ideals. 
  - `RightIdealClasses(OO)` in Magma
  - Cardinality is the class number





- Hecke does one split prime.  Can we use `PolyMake.jl` for totally positive and more general? 


## More to do
- Stabilizer/graph states 
- Parametric models / exponential families
- Solving pentagon and hexagon equations 
- Ross-Selinger and extensions
- Schur transform


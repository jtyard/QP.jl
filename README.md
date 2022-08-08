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
The invertible $R$-lattices in a central simple algebra over Frac$(R)$ and connecting 
$R$-orders form a groupoid.  

- 2-sided ideal class group
  - `TwoSidedIdealClassGroup(OO)` in Magma
- Set $\mathrm{Pic}_\ell(\mathcal{O})$ of left-equivalence classes of invertible right ideals. 
  - `RightIdealClasses(OO)` in Magma
  - Class number $h(\mathcal{O})$ is the number of right ideal classes
- Number of conjugacy classes, or types of orders
  - `#ConjugacyClasses(OO)` in Magma
- Hecke does one split prime.  Can we use `PolyMake.jl` for totally positive and more general? 


## More to do
- Stabilizer/graph states 
- Parametric models / exponential families
- Solving pentagon and hexagon equations 
- Ross-Selinger and extensions
- Schur transform


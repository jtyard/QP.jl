# QP.jl - Quantum Programming

- Abstract framework for modeling quantum systems
- Explicit computing in bases (vector, operator, etc)
- Implementation-agnostic (native julia, arb, FLINT, ..., TensorFlow)
- Mainly Julia/Oscar/Nemo/Hecke otherwise Magma, Sage, Pari

## Current capabilities

- Rings of polynomial functions on matrices over number fields 
- Useful matrices and tensor products 
- Qudit generalized Paulis
- Prime qudit generalized Clifford = Weil representation
- Constructing and investigating algebraic geometry of SIC-POVMs


## TODO 
- Generalize Weil to composite $N$
  -- Sort out even case
  -- Fix this [Oscar issue](https://github.com/oscar-system/Oscar.jl/issues/649) to implement $\mathrm{SL}(2,\mathbb{Z}/N)$ for all $N$.


## Towards arithmetic of quantum circuits
- 2-sided ideal class group
  - `TwoSidedIdealClassGroup(OO)` in Magma
- Set $\mathrm{Pic}_\ell$ of equivalence classes of invertible right ideals under left equivalence. 
  - `RightIdealClasses(OO)`
- Number of conjugacy classes, or types (Clark) of orders
  - `#ConjugacyClasses(OO)`
- Hecke does one split prime.  Can we use `PolyMake.jl` for totally positive and more general? 


## Longer term
- Stabilizer/graph states 
- Parametric models / exponential families
- Example: Solving pentagon and hexagon equations



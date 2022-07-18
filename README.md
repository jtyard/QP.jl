# QP.jl - Quantum Programming

Goals / desires
- Abstract framework for modeling quantum systems
- Explicit computing in bases (vector, operator, etc)
- Implementation-agnostic (native julia, arb, FLINT, ..., TensorFlow)
- Mainly Julia/Oscar/Nemo/Hecke otherwise Magma, Sage, Pari

## TODO 
- Generalize Weil to composite $N$
  -- Sort out even case
  -- Fix this [Oscar issue](https://github.com/oscar-system/Oscar.jl/issues/649) to implement $\mathrm{SL}(2,\mathbb{Z}/N)$ for all $N$.


## Towards arithmetic of quantum circuits
The following are already implemented in Hecke for not-totally-definite algebras but algorithms for definite algebras are lacking. 
- 2-sided ideal class group
  - `TwoSidedIdealClassGroup(OO)` in Magma
- Set $\mathrm{Pic}_\ell$ of equivalence classes of invertible right ideals under left equivalence. 
  - `RightIdealClasses(OO)`
- Number of conjugacy classes, or types (Clark) of orders
  - `#ConjugacyClasses(OO)`
- Hecke does one split prime.  `PolyMake.jl` for totally positive and more general? 

Explore this more


## Longer term
- Stabilizer/graph states 
- Parametric models / exponential families
- Example: Solving pentagon and hexagon equations



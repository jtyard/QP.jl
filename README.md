# QP.jl - for Quantum Programming

- 2-sided ideal class group
 - `TwoSidedIdealClassGroup(OO)` in magma
- Set $\mathrm{Pic}_\ell$ of equivalence classes of invertible right ideals under left equivalence. 
 - `RightIdealClasses(OO)`
- Number of conjugacy classes, or types (Clark) of orders
 - `#ConjugacyClasses(OO)`



## Goals / desires
- Abstract modeling of quantum systems
- Ability to define large, even infinite, tensor products and Hilbert spaces
- Work with whatever kinds of numbers you want (julia, arb, FLINT, ..., TensorFlow)
- Explicit computing in bases 
- Fast computation of stabilizer/graph states 
- Parametric models / exponential families


## Circuit synthesis over arithmetic groups
- Circuit synthesis for qubits with quaternion algebras and beyond
- Need to use `PolyMake.jl` for more general orders? 
- Solving pentagon and hexagon equations?


## Informal proof systems?

*Informal* proof system as workflow for searching for proofs using machine learning.  From MATH 239 to current research in physics and mathematics.
- How to represent information?
  - Scraping arXiv / hep-th
  - *concept classifier* 
  - Can a computer "discover" math or physics?
  - A possible answer: Computer very good at doing calculations.  In what sense could we possibly ask a computer to "understand" the result?  Can it learn to recognize / read / learn from existing proofs?
- Supervised / unsupervised learning: 
  - Supervised: Teach it
  - Unsupervised: Computers can easily generate examples by doing calculations.   
 Like, have it look at a paper and try to verify little parts as a check on its understanding.  
  - Publishing - How does it record its understanding or otherwise label what it finds "interesting".
  - In addition to julia also consider rust and go.








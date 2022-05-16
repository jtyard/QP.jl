# QP.jl - a package for quantum probability and quantum field theory

## Desired features
- Abstract modeling of quantum systems
- Ability to define large, even infinite, tensor products and Hilbert spaces
- Work with whatever kinds of numbers you want (julia, arb, FLINT, ..., TensorFlow)
- Explicit computing in bases 
- Fast computation of stabilizer/graph states 
- Parametric models / exponential families


### Circuit synthesis over arithmetic groups
- Circuit synthesis for qubits with quaternion algebras and beyond
- Is PolyMake.jl gonna help with the maximal orders?  How much is coded up by others???
- Solving pentagon and hexagon equations?
- A good first step is to code up examples from Vadym paper

### Tensor categories, $n$-categories, HoTT
- Proof theorists complain that julia doesn't have parameterized types in the strict sense of type theory

https://github.com/JuliaLang/julia/issues/6113

- There are however categorical logic systems built in julia.  One is 
https://github.com/epatters/Catlab.jl though I think there is at least one more.

- Monoidal categories are 2-categories.  Is there a set of "type definitions" for n-categories?  What about cobordisms or even just simplicial sets?

### Dependent types

https://homotopytypetheory.org/book/ (To do: solve Ch.9 exercises on categories and 2-categories so I can have some idea of what I'm talking about)

https://ncatlab.org/nlab/show/dependent+type+theory

https://ncatlab.org/nlab/show/type-theoretic+definition+of+category

### Homotopy Type Theory (HoTT)

Really should just solve the homework exercises in Chapter 9 of

HoTT is implemented in various automated proof systems

- Coq https://github.com/coq/coq is based on the Calculus of Inductive Constructions

- Agda https://wiki.portal.chalmers.se/agda/pmwiki.php and Idris  http://docs.idris-lang.org/en/latest/index.html are based on the Martin-Lof intensional theory of types.  Got Agda to install, seems popular enough, is built on Haskell (so I got that too now).  One more line to install Idris  ```cabal update; cabal install idris```,  which is supposed to be more of a full-blown programming language than Agda.  Presumably each can call Haskell code but I think this directly calls C.  However it seems less developed.

- https://arend-lang.github.io/ Is directly based on a cubical type theory https://ncatlab.org/nlab/show/cubical+type+theory which is a variant of HoTT in which the univalence axiom https://ncatlab.org/nlab/show/univalence+axiom is a theorem.

- https://leanprover.github.io/ (tried this breifly was annoyed that it uses git by default, taking over my vscode workspace)  Also seems HoTT is no longer supported and it seems sort of a waste of time to use it for this. 

### Martin-Lof dependent type theories
- Types = propositions and programs = proofs
- Deductive system: *rules* for deriving *judgements*
- judgement $a:A$ means "$a$ is an instance of the type $A$", or 
"$a$ is a witness, or evidence, for the truth of the proposition $A$". 
- Judgement $a:A$ of a proposition $A$ analogous to the set-theoretic proposition $a \in A$ but it is a judgement.  In actuality, $a:B$
- *Propositional equality* For $a,b:A$, there is an equality *proposition* $a=_Ab$
- *Definitional* or *judgemental* equality: There is also an equality judgement $a \equiv b :A$ 

### Homotopy types
- $A$ ia a *1-type* if for all $a,b:A$ all $p,q:x=y$ and $r,s:p=q$

### Univalence

### Other references

Displayed Categories https://arxiv.org/abs/1705.04296 (implemented in Coq)

### Informal proof systems?

Another related idea I'm having is that an *informal* proof system in julia could have advantages like the ability to use machine learning to search for proofs - in a sense it is teaching the computer to program, i.e. to find a way to transform one proposition into another. 

- How to get it to search for true things?  
- I really have practical things in mind, like 
  - MATH 239
  - calculus
  - physics
- Comes back to the bigger question of how to represent information?
  - Scraping arXiv 
  - *Concept classifier* (whatever that's gonna be)
  - Can a computer "discover" math or physics?
  - A possible answer: Computer very good at doing calculations.  In what sense could we possibly ask a computer to "understand" the result?  Can it learn from existing proofs?
- The idea: computer can augment its understanding by doing little calculations.  Like, have it look at a paper and try to verify little parts as a check on its understanding.  
  - Publishing - How does it record its understanding?
  - In long run don't get too attached to julia - also consider rust and go (especially if one of them has a better type system - what's the other one I found???)









# Type theories

## Tensor categories, $n$-categories, HoTT
- Proof theorists $n$ complain that julia doesn't have parameterized types in the strict sense of type theory

https://github.com/JuliaLang/julia/issues/6113

- There are however categorical logic systems built in julia.  One is 
https://github.com/epatters/Catlab.jl though I think there is at least one more.

- Monoidal categories are 2-categories.  Is there a set of "type definitions" for n-categories?  What about cobordisms or even just simplicial sets?

## Dependent types

https://homotopytypetheory.org/book/ (To do: solve Ch.9 exercises on categories and 2-categories so I can have some idea of what I'm talking about)

https://ncatlab.org/nlab/show/dependent+type+theory

https://ncatlab.org/nlab/show/type-theoretic+definition+of+category

## Homotopy Type Theory (HoTT)

Really should just solve the homework exercises in Chapter 9 of

HoTT is implemented in various automated proof systems

- Coq https://github.com/coq/coq is based on the Calculus of Inductive Constructions

- Agda https://wiki.portal.chalmers.se/agda/pmwiki.php and Idris  http://docs.idris-lang.org/en/latest/index.html are based on the Martin-Lof intensional theory of types.  Got Agda to install, seems popular enough, is built on Haskell (so I got that too now).  One more line to install Idris  ```cabal update; cabal install idris```,  which is supposed to be more of a full-blown programming language than Agda.  Presumably each can call Haskell code but I think this directly calls C.  However it seems less developed.

- https://arend-lang.github.io/ Is directly based on a cubical type theory https://ncatlab.org/nlab/show/cubical+type+theory which is a variant of HoTT in which the univalence axiom https://ncatlab.org/nlab/show/univalence+axiom is a theorem.

- https://leanprover.github.io/ (tried this breifly was annoyed that it uses git by default, taking over my vscode workspace)  Also seems HoTT is no longer supported and it seems sort of a waste of time to use it for this. 

## Univalence

## Other references

Displayed Categories https://arxiv.org/abs/1705.04296 (implemented in Coq)
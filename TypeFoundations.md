# Type foundations

## Formal language

A [formal language](https://en.wikipedia.org/wiki/Formal_language) is a subset $L \subset \Sigma^*$ of strings over a (usually finite) alphabet $\Sigma$.

Kleene hierarchy
- $\mathsf{ALL}$ all languages
- $\mathsf{CS}$ context-sensitive languages
- $\mathsf{CF}$ context-free languages
- $\mathsf{RE}$ recursively enumerable languages

## Types
- [Type theories](https://ncatlab.org/nlab/show/type+theory) are formal languages consisting of terms.  Some of these terms are called types and each term has a definite type.  So is it a formal language $L$ together with a map 
$$\mathsf{type}:L \to L?$$
Applications include
    - Syntax for category theory
    - Formal verification and automated theorem proving 
    - Foundations of mathematics
    - Logic, proof theory, computability theory
- [Propositions as types](https://ncatlab.org/nlab/show/propositions+as+types): Regard each type $A$ as a proposition each term $t:A$ as a proof, or witness, of $A$.    


## Extensional vs intensional
- [Extensional](https://ncatlab.org/nlab/show/extensional+type+theory) [vs](https://ncatlab.org/nlab/show/type+theory#ExtensionalIntensional) [Intensional](https://ncatlab.org/nlab/show/intensional+type+theory)
## Dependent types
- [dependent type](https://ncatlab.org/nlab/show/dependent+type) 
- [dependent type theory](https://ncatlab.org/nlab/show/dependent+type+theory)
- A morphism $f : A \to B$ in a category gives a dependent type $f(a):B$ 
- Can [define] categories (https://ncatlab.org/nlab/show/type-theoretic+definition+of+category) using dependent type theories.  Similarly with 2-categories, of which monoidal categories are an example.  What is known about type definitiosn of $n$-categories? What about cobordisms or even just simplicial sets?
    -To do: solve Ch.9 exercises in [HoTT](https://homotopytypetheory.org/book/) on categories and 2-categories so I can have some idea of what I'm talking about.




## Martin-Lof dependent type theory 
- [Martin-Lof dependent type theory](https://ncatlab.org/nlab/show/Martin-L%C3%B6f+dependent+type+theory)


## Homotopy Type Theory (HoTT)

### Univalence
- [Univalence axiom](https://ncatlab.org/nlab/show/univalence+axiom)

### Implementations of HoTT
HoTT is implemented in various automated proof systems

- Coq https://github.com/coq/coq is based on the Calculus of Inductive Constructions

- Agda https://wiki.portal.chalmers.se/agda/pmwiki.php and Idris  http://docs.idris-lang.org/en/latest/index.html are based on the Martin-Lof intensional theory of types.  Got Agda to install, seems popular enough, is built on Haskell (so I got that too now).  One more line to install Idris  ```cabal update; cabal install idris```,  which is supposed to be more of a full-blown programming language than Agda.  Presumably each can call Haskell code but I think this directly calls C.  However it seems less developed.

- https://arend-lang.github.io/ Is directly based on a cubical type theory https://ncatlab.org/nlab/show/cubical+type+theory which is a variant of HoTT in which the univalence axiom https://ncatlab.org/nlab/show/univalence+axiom is a theorem.

- https://leanprover.github.io/ (tried this breifly was annoyed that it uses git by default, taking over my vscode workspace)  Also seems HoTT is no longer supported and it seems sort of a waste of time to use it for this. 
## Julia
- Proof theorists complain that julia doesn't have parameterized types in the strict sense of type theory

https://github.com/JuliaLang/julia/issues/6113

- There are however categorical logic systems built in julia.  One is 
https://github.com/epatters/Catlab.jl though I think there is at least one more.


## Other references

Displayed Categories https://arxiv.org/abs/1705.04296 (implemented in Coq)
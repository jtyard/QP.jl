# Type foundations

## Formal language

- A [formal language](https://en.wikipedia.org/wiki/Formal_language) is a subset $L \subset \Sigma^*$ of strings over a (usually finite) alphabet $\Sigma$.
- Formal system is a formal language along with... 
- Formal logic specifies a subset of $L$ to be used as [judgements](https://ncatlab.org/nlab/show/judgment)

Kleene hierarchy
- $\mathsf{ALL}$ all languages
- $\mathsf{CS}$ context-sensitive languages
- $\mathsf{CF}$ context-free languages
- $\mathsf{RE}$ recursively enumerable languages

## Types
- [Type theories](https://ncatlab.org/nlab/show/type+theory) are formal languages consisting of terms.  Some of these terms are called types and each term has a definite type.  So is it a formal language $L$ together with a map 
$$\mathsf{type}:L \to L?$$
- In set theory we use the propositional calculus to build sets out of sets, whereas type theory has only types i.e. we use types to build types from types.  
- Judgements belong to the meta-language.
- Write $t:A$ for the judgment "$t$ is a term of type $A$". 
- [Propositions as types](https://ncatlab.org/nlab/show/propositions+as+types): Identify each type $A$ with the proposition "there is a term $t:A$".  Term $t$ of type $A$ viewed as a proof, or witness, of $A$ (or rather, that there is a term $t:A$.)
    - Generalizes propositional logic, which can be viewed as rules for doing things like turning a proof of $A$ and a proof of $B$ into a proof of $A \wedge B$.
- Applications include
    - Syntax for category theory
    - Formal verification
    - Automated theorem proving 
    - Foundations of mathematics
    - Logic, proof theory, computability theory
 


## Extensional vs intensional
- [Extensional](https://ncatlab.org/nlab/show/extensional+type+theory) [vs](https://ncatlab.org/nlab/show/type+theory#ExtensionalIntensional) [Intensional](https://ncatlab.org/nlab/show/intensional+type+theory)
## Dependent types
- [dependent type](https://ncatlab.org/nlab/show/dependent+type) 
- [dependent type theory](https://ncatlab.org/nlab/show/dependent+type+theory)
- A morphism $f : A \to B$ in a category gives a dependent type $f(a):B$ 
- Can [define] categories (https://ncatlab.org/nlab/show/type-theoretic+definition+of+category) using dependent type theories.  Similarly with 2-categories, of which monoidal categories are an example.  What is known about type definitiosn of $n$-categories? What about cobordisms or even just simplicial sets?
    -To do: solve Ch.9 exercises in [HoTT](https://homotopytypetheory.org/book/) on categories and 2-categories so I can have some idea of what I'm talking about.




## Martin-Lof type theory 
- [Martin-Lof type theory](https://en.wikipedia.org/wiki/Intuitionistic_type_theory) aka intuitionistic type theory or constructive type theory
- [Martin-Lof dependent type theory](https://ncatlab.org/nlab/show/Martin-L%C3%B6f+dependent+type+theory) aka 

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


## Chapter -1?

A *language* $L$ over an *alphabet* $A$ is just a subset $L \subset A^*$.

A *formal language* is a language accepted by some regular expression / formal grammar / finite automaton. 

context free $\subset$ context sensitive $\subset$ recursively enumerable $\subset$ computable $\subset \mathsf{ALL} = 2^{\\{0,1\\}^*}$

Maybe some complexity classes?

## Martin-L\"{o}f dependent type theories
A deductive system, such as first-order predicate logic, applies 
- Types = propositions and programs = proofs

- Deductive system: *rules* for deriving *judgements*

- judgement $a:A$ means `'$a$ is an instance of the type $A$'', or 
`'$a$ is a witness, or evidence, for the truth of the proposition $A$''. 

- Judgement $a:A$ of a proposition $A$ analogous to the set-theoretic proposition $a \in A$ but it is a judgement.  In actuality, $a:B$

- *Propositional equality* For $a,b:A$, there is an equality \emph{proposition} $a=_Ab$

- *Definitional* or *judgemental equality: There is also an equality judgement $a \equiv b:A$ 

*Homotopy types*: $A$ ia a \emph{1-type} if for all $a,b:A$ all $p,q:x=y$ and $r,s:p=q$



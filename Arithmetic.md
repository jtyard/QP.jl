# Arithmetic of quantum circuits
Let $R$ be an integral domain and $B$ a central simple algebra over $\mathrm{Frac}(R)$.  
- An **$R$-lattice** is a finitely generated torsion-free $R$-module.
- A **fractional $R$-ideal** $I$ of $B$ is a full-rank $R$-sublattice.
  - **Left order** $`O_\ell(I) = \{ b \in B : bI \subset I \}`$
  - **Right order**  $`O_r(I) = \{ b \in B : bI \subset I \}`$
  - $I$ is an $O_\ell(I)-O_r(I)$-bimodule, a left fractional $O_\ell(I)$ 
ideal and a right fractional $O_r(I)$-ideal.
  - There is a **category** whose objects are the orders and whose morphisms are the fractional ideals, with composition given by proper products (requiring matching left and right orders).  Fractional ideals **connect** their left/right orders.
- $O\simeq O'$ iff $bOb^{-1} = O'$ for some $b \in B^\times$ 
iff they are connected by a principal ideal $bO = O'b$. So the isomorphism classes are the conjugacy classes.
- An **invertible ideal** of $B$ is a fractional ideal $I$ with an **inverse** $I^{-1}$, satisfying 
$I^{-1}I = O_R(I)$ and therefore $I I^{-1} = O_\ell(I)$.
  - The orders and the invertible ideals form a groupoid called the **core** of the above category.
- A **normal ideal** of $B$ is a fractional $R$-ideal with maximal left/right orders.
  - normal $\Rightarrow$ invertible.  
  -  The maximal orders and normal ideals of $B$ form the **Brandt groupoid**.
- An **integral ideal** of $B$ is a normal $R$-ideal $I$ with $I \subset O_\ell(I) \cap O_r(I)$. 
- A **maximal left ideal** $I$ of $B$ is a maximal left ideal of $O_\ell(I)$.  
- Each integral ideal $I$ of $B$ can be factored into a **proper product** $I = M_1 M_2 \cdots M_n$ of maximal left ideals $M_i$ of $B$ with $O_r(M_i) = O_\ell(M_{i+1})$ for $1 \leq i < n$, 
where $n$ is the length of $O_\ell(I)/I$.


## Wish list

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
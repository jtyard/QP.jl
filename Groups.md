# Groups in Oscar/GAP







## Group actions in Oscar

Group action definitions and examples
- Oscar/src/Groups/action.jl
- Oscar/test/Groups/action.jl

What are gsets for vs just defining the action as functions?
- Oscar/src/Groups/gsets.jl
- Oscar/test/Groups/gsets.jl


There are many ways to act on polynomials with groups.  Mainly of interest is acting on Nemo.MPolyElem with PermGroupElem or MatrixGroupElem.

These are defined in Oscar/src/Groups/action.jl

  on_indeterminates(f::GAP.GapObj, p::PermGroupElem)
  on_indeterminates(f::Nemo.MPolyElem, p::PermGroupElem)
  on_indeterminates(f::MPolyIdeal, p::PermGroupElem)
  on_indeterminates(f::GAP.GapObj, p::MatrixGroupElem)
  on_indeterminates(f::Nemo.MPolyElem{T}, p::MatrixGroupElem{T, S}) where T where S
  on_indeterminates(f::MPolyIdeal, p::MatrixGroupElem)

## Representations

Given a representation G -> GL(n,F), get action G -> R[G]


## 




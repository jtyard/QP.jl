using QP 
using Oscar


# Test cases for groups 

is_isomorphic(SL(2,2),symmetric_group(3))

is_isomorphic(PSL(2,3),alternating_group(4))

is_isomorphic(PSL(2,4),symmetric_group(4))

is_isomorphic(PSL(2,5),alternating_group(5))

is_isomorphic(PGL(2,5),symmetric_group(5))

[dim(V) for V in irreducible_modules(SL(2,ZN(4)))]
###############
# Polynomials in matrix entries with action of (Z/N)^2 ⋊ SL_2(Z/2N) and related groups
# So far nothing good here just playing with types
###############

using Oscar
using Memoize 

# Eventually implement types like this
abstract struct MatPolyRing <: MPolyRing
abstract struct MatPoly <: MPolyElem




    
# Coordinates for operator basis of matrix elements Q[X], X = sum_ij Xij Eij  

# Coordinates Q[Y] = sum_ij Y_ij Delta_ij 

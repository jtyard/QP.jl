using Oscar

# Miller-Sturmfels - Combinatorial commutative Algebra

# We create the graded ring from example 8.2 
A = abelian_group([0,2])
S, _ = PolynomialRing(QQ,3)
SA, (x1,x2,x3) = grade(S,[A([1,1]),A([-2,1]),A([1,0])])
# The following doesn't work
#homogeneous_component(R2A,A([2,0]))






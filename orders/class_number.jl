# Class number and related methods for quaternion orders

using Oscar, QP  


K,s = quadratic_field(2)
OK = maximal_order(K)
A = quaternion_algebra(QQ,-1,-1)
B = quaternion_algebra(K,-1,7(1+s))

OA = maximal_order(A)
OB = maximal_order(B)

a,b = B.std

discriminant(B)
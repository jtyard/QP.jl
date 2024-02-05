using Oscar
using QP

export quaternion_algebra

# A wrapper that coerces the generators automatically - should merge into OSCAR
quaternion_algebra(K,a,b) = QuaternionAlgebra(K,K(a),K(b))

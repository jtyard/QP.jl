# SICs on Calabi-Yaus
# Verifying that each SIC lies on a unique X_1^N + ... + X_N^N - e^{t/N} X_1 ... X_N
# with e^{t/N} a unit in the complex ring class field of the corresponding quadratic order.
# 
# In particular this should be real when there are real fiducials.

using Oscar, QP

function cyparam(v)
    N = length(v)
    sum([a^N for a in v])//prod([a for a in v])
end

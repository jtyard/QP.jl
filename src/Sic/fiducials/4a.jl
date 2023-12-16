S4 = SicData(4,build_nf=true)
F = S4.F
#F = number_field(S4.rcf,using_stark_units = true)

# Define some constants
w2 = sqrt(F(2))
w5 = sqrt(F(5))
w10 = w2*w5
r1=sqrt(w5+1)
I=sqrt(F(-1))
phi=S4.F[
    8
    ((w10+w2-2*w5-2)*r1-4)*I+(w10+w2)*r1+4*w2-4
    (8*w2-8)*I
((w10+w2-2*w5-2)*r1+4)*I+(w10+w2)*r1-4*w2+4]

#gal, sig = automorphism_group(S4.rcf)
#abs2c(x) = x*c(x)

# Rank-1 density matrix 
Phi = phi*dagger(phi,S4.c); 
return trace(Phi)^-1 * Phi
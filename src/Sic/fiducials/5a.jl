S5 = SicData(5, build_nf=true)
F = S5.F

# Define some constants
w3 = sqrt(F(3));
w5 = sqrt(F(5));
w15 = w3*w5;
I = sqrt(F(-1));

r1=sqrt(1//2*w5+5//2)
r2=sqrt((-1//8*w3+1//8*w5+3//8)*r1+1//8*w3*w5+5//8*w3+1//16*w5-5//16)

# Find a complex conjugation 
#gal, sig = automorphism_group(S5.rcf.A)
#possible_c = findall([map(sig(g),[I,w3,w5,r1,r2]) == [-I,w3,w5,r1,r2] for g in gal])
#println("Better only be one choice in this list " * string(possible_c))
#c = sig(gal[possible_c[1]])
#gal, sig = automorphism_group(S5.rcf)
#c = complex_conjugation(S5.rcf, S5.inf[2])
#abs2c(x) = x*c(x)

# The fiducial vector
phi = F[16*r1
    (((8*w3-10*w5-2)*r1+20*w3-16*w5)*r2+((2*w5+2)*r1+(2*w5+10)))*I+((-6*w15-6*w3+8*w5+24)*r1+(-12*w15-20*w3+12*w5+40))*r2+(4*w5-8)*r1-4*w15-4*w5
    (((16*w3-14*w5+6)*r1+(-2*w15+30*w3-24*w5))*r2+((w15-5*w3+4)*r1-3*w5+5))*I+((2*w15-2*w3+4*w5-12)*r1-2*w5-10)*r2+(w5+3)*r1+3*w15-5*w3-2*w5+10
    (((-6*w15-6*w3+6*w5+14)*r1+(-8*w15-20*w3+16*w5+20))*r2+((w15+5*w3-2*w5+6)*r1+4*w15-2*w5))*I+((2*w15+2*w3-2*w5-18)*r1+20*w3-12*w5)*r2+(3*w5-1)*r1-2*w15
    (((-8*w15+12*w3-12*w5+32)*r1+(-10*w15+10*w3-6*w5+30))*r2+((-2*w15+w5+3)*r1+w15+5*w3+3*w5+5))*I+((-4*w3-4)*r1+2*w15-10*w3+2*w5-10)*r2+(-w15+5*w3+2*w5)*r1+3*w15+5*w3+w5+5] 

# Rank-1 density matrix 
Phi = phi*dagger(phi,S5.c); 
Phi = trace(Phi)^-1 * Phi
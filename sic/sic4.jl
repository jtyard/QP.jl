using QP
using Oscar 

N = 4

S4 = SicData(N,build_nf=true)

K = S4.K 
F = S4.rcf.A

w2= -sqrt(F(2))
w5= -sqrt(F(5))
w10=w2*w5
r1=sqrt(w5+1)
I=sqrt(F(-1))

phi=F[
8
((w10+w2-2*w5-2)*r1-4)*I+(w10+w2)*r1+4*w2-4
(8*w2-8)*I
((w10+w2-2*w5-2)*r1+4)*I+(w10+w2)*r1-4*w2+4
]


# Find a complex conjugation 
gal, sig = automorphism_group(S4.rcf)
c = complex_conjugation(S4.rcf, S4.inf[1])
#possible_c = findall([map(sig(g),[I,w2,w5,r1]) == [-I,w2,w5,r1] for g in gal])
#println("Better only be one choice in this list " * string(possible_c))
#c = sig(gal[possible_c[1]])
abs2c(x) = x*c(x)

# Rank-1 density matrix 
Phi = phi*transpose(map(c,phi)); Phi = trace(Phi)^-1 * Phi

println(is_fiducial(Phi))

#phi = Matrix(d,1,philist);

#C,_ = cyclotomic_field(8)
#C_to_F = hom(C,F,zetaN(8,F))
#
## compute the overlaps
#overlaps = [[trace(map(C_to_F,heis(Z4[i j]))*Phi) for j in 0:3] for i in 0:3]
#
##and check that they all equal 1/(d+1)
#map2(f,a) = map(x->map(f,x),a) # so map2(f,[[ ]]) = [[ f(..) ]]
#display(map2(abs2c,overlaps))
#
##Check that the fiducial satisfies the harmonic polyomial equations directly
#println()
#display([[evaluate(h(Z4[i j]),vec(Phi)) for j in 0:3] for i in 0:3])



R = QQX(N,graded=true); 
Rp = QQXp(N,graded=true); 
I = Im(N,graded=true) + Ih(N,graded=true) + Ic(N,graded=true)
#Rt = QQXt(N); 
S = ProjectiveScheme(R,I);

#pd = primary_decomposition(I)
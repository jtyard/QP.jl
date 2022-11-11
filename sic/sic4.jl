using QP
using Oscar 

#Z4 = ZN(4)
S4 = SicData(4,build_nf=true)

K = S4.K 
F = S4.rcf.A

w2= -sqrt(F(2))
w5= -sqrt(F(5))
w10=w2*w5
r1=sqrt(w5+1)
I=sqrt(F(-1))

phiv=F[
8
((w10+w2-2*w5-2)*r1-4)*I+(w10+w2)*r1+4*w2-4
(8*w2-8)*I
((w10+w2-2*w5-2)*r1+4)*I+(w10+w2)*r1-4*w2+4
]


# Find a complex conjugation 
gal, sig = automorphism_group(S4.rcf.A)
possible_c = findall([map(sig(g),[I,w2,w5,r1]) == [-I,w2,w5,r1] for g in gal])
println("Better only be one choice in this list " * string(possible_c))
c = sig(gal[possible_c[1]])
abs2c(x) = x*c(x)


#F = number_field(S.rcf,using_stark_units = true)
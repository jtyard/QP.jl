using QP 
using Oscar

S5 = SicData(5,build_nf=true)
F = S5.rcf.A

# Define some constants
w3 = sqrt(F(3));
w5 = sqrt(F(5));
w15 = w3*w5;
I = sqrt(F(-1));

r1=sqrt(1//2*w5+5//2)
r2=sqrt((-1//8*w3+1//8*w5+3//8)*r1+1//8*w3*w5+5//8*w3+1//16*w5-5//16)

# Find a complex conjugation 
gal, sig = automorphism_group(S5.rcf.A)
possible_c = findall([map(sig(g),[I,w3,w5,r1,r2]) == [-I,w3,w5,r1,r2] for g in gal])
println("Better only be one choice in this list " * string(possible_c))
c = sig(gal[possible_c[1]])


set_verbose_level(:ClassField, 2)


# end(codomain(f)) -> end(domain(f)) 
function map_end(g,f)
    return hom(domain(f),domain(f),[inv(f)(g(f(a))) for a in gens(domain(f))]...)
end


Fs, Fs_to_F = absolute_simple_field(F)
Fss, Fss_to_Fs = simplify(Fs)
Fss_to_F = compose(Fss_to_Fs,Fs_to_F)

css = map_end(c,Fss_to_F)

C60,z60 = cyclotomic_field(60)
z3 = z60^20
z4 = z60^15
z5 = z60^12
C60_to_Fss = hom(C60,Fss,zetaN(60,Fss))
Frel, Frel_to_Fss = relative_simple_extension(Fss,C60)
e = discriminant(Frel) 
Cgal, Csig = automorphism_group(C60)
possible_sig2 = findall([map(Csig(g),[z3,z4,z5]) == [z3,z4,z5^2] for g in Cgal])
possible_sig3 = findall([map(Csig(g),[z3,z4,z5]) == [z3,z4,z5^3] for g in Cgal])
possible_sig4 = findall([map(Csig(g),[z3,z4,z5]) == [z3,z4,z5^4] for g in Cgal])
possible_tau3 = findall([map(Csig(g),[z3,z4,z5]) == [z3^-1,z4,z5] for g in Cgal])
possible_tau4 = findall([map(Csig(g),[z3,z4,z5]) == [z3,z4^-1,z5] for g in Cgal])
sig2 = Csig(Cgal[possible_sig2[1]])
sig3 = Csig(Cgal[possible_sig3[1]])
sig4 = Csig(Cgal[possible_sig4[1]])
tau3 = Csig(Cgal[possible_tau3[1]])
tau4 = Csig(Cgal[possible_tau4[1]])

e*sig4(e)*sig2(e)*sig3(e)
# After all this I can finally say that F = Q(z60,sqrt(e)) for e = -z_60^14 + z_60^10 + z_60^8 - z_60^6 - 2*z_60^4 + 1
# Note that norm_Q(z12)(e) = 5.
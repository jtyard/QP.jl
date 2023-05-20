#

using QP
using Oscar

N = 3; 

gr = true
R = QQX(N,graded=gr); 
Rp = QQXp(N,graded=gr); 
I = Im(N,graded=gr) + Ih(N,graded=gr); 
#Rt = QQXt(N); 

Isat = saturation(Rp,I)

S = ProjectiveScheme(R,I);
#println("Computing primary decomposition")


pd = primary_decomposition(I,alg=:SY)
P = [p[2] for p in pd]
#Rs, ph = quo(R,Im(N)); 

#F = cyclotomic_field(12)[1]

#rational_solutions(change_base_ring(F,I)) 


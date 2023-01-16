#
using QP
using Oscar

N = 2; 

R = QQX(N,graded=true); 
Rp = QQXp(N,graded=true); 
I = Im(N,graded=true) + Ih(N,graded=true); 
#Rt = QQXt(N); 
S = ProjectiveScheme(R,I);

pd = primary_decomposition(I)

#Rs, ph = quo(R,Im(N)); 

F = cyclotomic_field(12)[1]

rational_solutions(change_base_ring(F,I)) 


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

# X = heis(0,1,2)
# Y = heis(1,1,2)
# Z = heis(1,0,2)


Phi = (identity_matrix(F,2) + map(F,heis(0,1,2) + heis(1,0,2) + heis(1,1,2))/sqrt(F(3)))/2

sic = heis_orbit(Phi)
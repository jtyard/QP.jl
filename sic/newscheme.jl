#
using Oscar, QP

N = 7

R = matrix_polynomial_ring(QQ,N)

Ih = Ihplus(R)
Im = Iminors(R,2)
It = Itorus(R,Z7[2 0; 0 4])
Ic = QP.Ic(R)

Iab = Ih + Im + It
Ib = Iab + Ic
#Sab = proj(Iab)
#Sb = proj(Ib)

#degree(Sab), degree(Sb)
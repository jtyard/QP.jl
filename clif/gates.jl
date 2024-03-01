# Some qubit gates
C3, z3 = cyclotomic_field(3)
C4, i = cyclotomic_field(4,"i")
C8, z8 = cyclotomic_field(8)
C9, z9 = cyclotomic_field(9)

C12,z12 = cyclotomic_field(12)

s3 = sqrt(C12(-3))
#z3 = z12^4
#i = z12^3


t12 = C4[3 1-i; 1+i -3] # Unnormalized super golden gate for PU(2)

X = QQ[0 1; 1 0]
Y = C4[0 -i;i 0]
Z = QQ[1 0; 0 -1]

u = (1//2)*(1 + i*(X+Y+Z))


# Isometry onto an embedded qubit
U = C4[1+i 0; 0 1; 0 1]
(1//2)*U*dagger(U)



# WHy isn't this working???
[det((1//2)*dagger(U)*weil_U(g)*U) for g in SL(2,Z3)]

#Matrix entry of the fiducial (0 : 1 : -1)
r  = QQ[0; 1; -1]
[((1//2)*transpose(r)*weil_U(g)*r)[1] for g in SL(2,Z3)]
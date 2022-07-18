# 
N = 5; C12, z12 = cyclotomic_field(12); I = z12^3; z3 = z12^4; w3 = 1+2*z3; s3 = I*w3;

#s3 := Sqrt(F! 3);
#I := Sqrt(F!-1);
#w3 := I*s3;     

#Kinf := AbsoluteField(sub<F|w3,s3>);

ZC12 = maximal_order(C12); GKinf = galois_group(C12);

g2 = (1125//2)*(3+4*I)*(1-w3); g3 = 4125*s3*(11-2*I);

a = -g2//4; b = -g3//4;

E = EllipticCurve([a, b]);

"Computing minimal model";
Emin := MinimalModel(E);
jInvariant(E);
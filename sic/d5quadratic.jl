# Trying to realize the d=5 solution field as a "natural" quadratic extension of Q(zeta_60)

using Oscar

# Rationals as a number field
_, x = PolynomialRing(QQ)
Q,_ = NumberField(x-1,"a")
Z = maximal_order(Q)
RC = absolute_simple_field(NumberField(ray_class_field(60*Z)))[1]

K,s = quadratic_field(3)
OK = maximal_order(K)
FK = NumberField(ray_class_field(5*OK,[infinite_places(K)[1]]))
F,FK_to_F = absolute_simple_field(FK)

subs = subfields(F)
rc = [f for f in subs if degree(f[1]) == 8][1][1]
_,rc_to_RC = isisomorphic(rc,RC)

Frc,F_to_Frc = relative_simple_extension(F,rc)

dd = discriminant(Frc)

_,y = PolynomialRing(K,"y")
f = 1296y^8 - (648 - 1080s)y^7 + 648y^6 + (864 - 360s)y^5 - (540-360s)y^4 + (144-60s)y^3 + 18y^2 - (3-5s)y + 1



C60,z60 = CyclotomicField(60)

factor(ZZ(discriminant(C60)))



degree(ray_class_field(5*OK,infinite_places(K)))
factor(ZZ(absolute_discriminant(ray_class_field(5*OK,infinite_places(K)))))

dd = discriminant(ray_class_field(5*OK,infinite_places(K)))
factor(37252902984619140625)
factor(6103515625)

println(collect(factor(dd)))





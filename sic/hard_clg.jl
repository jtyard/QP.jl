using Oscar

_,x = PolynomialRing(QQ)
F = NumberField(x^32 + 5130446197409069928117407*x^16 + 6218559606652239161766175410604256580325035025921)[1]

set_verbose_level(:ClassGroup,3)
#set_verbose_level(:ClassGroup_time,1)

@time class_group(F)

#=
Indeed this is just a cyclotomic extension of the quadratic base field

julia> S = SicData(68)
SicData(68, 4485, 4485, 1, Number field over Rational Field with defining polynomial x^2 - 4485, Maximal order of Number field over Rational Field with defining polynomial x^2 - 4485 
with basis nf_elem[1, 1//2*s + 1//2], 1//2*s + 67//2, Order of Number field over Rational Field with defining polynomial x^2 - 4485
with Z-basis NfOrdElem[1, 1//2*s + 4485//2], Order of Number field over Rational Field with defining polynomial x^2 - 4485
with Z-basis NfOrdElem[1, 1//2*s + 2242], 1//2*s + 2242, Class field defined mod (<136, 136>, InfPlc[Real place of 
Number field over Rational Field with defining polynomial x^2 - 4485
corresponding to the root [-66.9701426010128785361058352970499431197457436493339321877271723905415953982421772224171794268016751025295220332918245837163799131333565067247000428859541 +/- 3.38e-152], Real place of 
Number field over Rational Field with defining polynomial x^2 - 4485
corresponding to the root [66.9701426010128785361058352970499431197457436493339321877271723905415953982421772224171794268016751025295220332918245837163799131333565067247000428859541 +/- 3.38e-152]]) of structure Abelian group with structure: (Z/2)^4 x Z/4 x Z/288)

julia> S.K
Number field over Rational Field with defining polynomial x^2 - 4485

julia> cyclotomic_extension(S.K,32)
Cyclotomic extension by zeta_32 of degree 32

julia> absolute_simple_field(cyclotomic_extension(S.K,32))
Number field over Rational Field with defining polynomial x^32 + 5130446197409069928117407*x^16 + 6218559606652239161766175410604256580325035025921

julia> class_group(S.K)
(GrpAb: (Z/2)^2, ClassGroup map of 
Set of ideals of Maximal order of Number field over Rational Field with defining polynomial x^2 - 4485 
with basis nf_elem[1, 1//2*s + 1//2]
)

julia> fundamental_discriminant(4485)
4485

julia> factor(4485)
1 * 5 * 13 * 3 * 23
=#
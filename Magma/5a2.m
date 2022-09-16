// 

als := {@ f2(m) : m in M | m ne 0 @};

al := als[1];

Fabs := AbsoluteField(F);

GalFQ, _ , aut := AutomorphismGroup(Fabs);

GalFK := sub<GalFQ | [g : g in GalFQ | aut(g)(w3) eq w3]>;

GalFKinf := sub<GalFK | [g : g in GalFK | aut(g)(I) eq I]>;

// ccgal := [g : g in GalFK | aut(g)(al) eq al][1];

cocyc := func<j,k|zd^(Integers() ! (3*j[2]*k[1] - 3*j[1]*k[2])) >;
pairing := func<j,k|zd^(Integers() ! (j[2]*k[1] - j[1]*k[2])) >;

// Need to fix this as I don't actually have a strongly centered fiducial.
// But note another way to prase it: they are ALL Stark units, up to roots of unity. 
a := al;

cc(a) eq a;

// Given automorphism of F/K, returns corresponding automorphism of Fabs/Q
absaut := func<g | hom<Fabs -> Fabs | x :-> Fabs ! ( g(F ! x))>>;

GalFLplus := sub<GalFK | [g : g in GalFK | aut(g)(a) eq a]>; 

GalFKplus := sub<GalFK | [g : g in GalFK | aut(g)(w3) eq w3 and aut(g)(I) eq I]>;

// Weird.  Seems Magma doesn't always embed the field the same way
// Sometimes need to change al by a root of unity to strongly center it.

Lplus := FixedField(Fabs,GalFLplus);
Kplus := FixedField(Fabs,GalFKplus);

VK, FtoVK := VectorSpace(Fabs,K);
VKplus, FtoVKplus := VectorSpace(Fabs,Kplus);
VK2, LplustoVK2 := VectorSpace(Lplus,K);

L1 := NumberFieldLattice([FtoVKplus(Fabs ! a) : a in als]);

L2 := NumberFieldLattice([LplustoVK2(Lplus ! Fabs ! a) : a in als]);




d:=5;
K := QuadraticField((d-3)*(d+1));
OK := RingOfIntegers(K);

"Creating ray class group";
rayclassgroup, representing_ideal := RayClassGroup(d*OK,[1,2]);

"Building abelian extension";
FA := AbelianExtension(representing_ideal);

"Constructing number field";
F:= NumberField(FA);

P2 := Factorization(2*OK)[1][1];
P3 := Factorization(3*OK)[1][1];

s3 := ArtinMap(FA)(P3);
s2 := ArtinMap(FA)(P2);

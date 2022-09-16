d:=4;
dd:=8;

K := QuadraticField((d-3)*(d+1));
OK := RingOfIntegers(K);

"Creating ray class group";
rayclassgroup, representing_ideal := RayClassGroup(dd*OK,[1,2]);

"Building abelian extension";
FA := AbelianExtension(representing_ideal);

"Constructing number field";
F:= NumberField(FA);

gal := func<s | ArtinMap(FA)(representing_ideal(s))>;

//Gal := [gal(s) : s in rayclassgroup];
//c := Gal[3];

// Find (or define) a primitive dth roots of unity
//"Computing ddth root of unity";
//zdd := F ! RootOfUnity(2*d,OptimizedRepresentation(AbsoluteField(F)));

//zd := zdd^2;

w2:= -Sqrt(F!2);
w5:= -Sqrt(F!5);
w10:=w2*w5;
r1:=Sqrt(w5+1);
I:=Sqrt(F!(-1));

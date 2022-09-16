d:=8;
K := QuadraticField((d-3)*(d+1));
OK := RingOfIntegers(K);
R := quo<OK|d*OK>;
randK := [Random(OK,100) : i in [1..10000] ];
quos := {<R!a,RealSigns(a)> : a in randK};
#quos;
4*d^2;


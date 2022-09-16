// A conjecture Q(sqrt(D_0)) contains an element of norm d/D_0

for d in [4..103] do
    D0 := FundamentalDiscriminant((d-3)*(d+1));
    K := QuadraticField(D0);
    if NormEquation(K,d/D0) then
	print d, D0, NarrowClassNumber(K), ClassNumber(K);
    end if;

end for;

D0 := FundamentalDiscriminant(5);
K:= QuadraticField(D0);
s5:= K.1;
u := (3+s5)/2;
for r in [1..5] do
    d := u^r + u^-r + 1;
    if NormEquation(K,d/D0) then
        print d, D0, ClassNumber(K);
    end if;

end for;


K := QuadraticField(2);
OK := RingOfIntegers(K);
s2 := K.1;

P:= (3+s2)*OK;
PP := (3-s2)*OK;

P*PP;

uf := 1+s2;
quo<OK|P> ! uf;
quo<OK|PP> ! uf;
quo<OK|P> ! s2;
quo<OK|PP> !s2;
//SL2OK := SL(2,OK);

SL2 := SL(2,7);
GL2 := GL(2,7);

G := SL2;
T := sub<G | G ! [5,0,0,3]>;
B := sub<G | T, G ! [1,1,0,1]>;

for B1 in Conjugates(G,B) do
    for B2 in Conjugates(G,B) do
	S := B1 meet B2;
	if not IsIsomorphic(S,T) then
	    print #S;
	end if;
    end for;
end for;

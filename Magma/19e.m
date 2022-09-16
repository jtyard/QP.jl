d  := 19
K  := QuadraticField(d);
OK := RingOfIntegers(K);
s5 := K.1;

P  := (7+8*s5)*OK;
PP := (7-8*s5)*OK;

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

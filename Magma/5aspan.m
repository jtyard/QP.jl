phis := {};
for m in M do
    phis := phis join {Delta(m)*phi};
end for;

V := VectorSpace(F,d);

Min([Rank(sub<V | [V ! Eltseq(phi) : phi in subs]  >) : subs in Subsets(phis,5)]);

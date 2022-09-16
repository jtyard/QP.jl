
phis := {};
for m in M do
    phis := phis join {Delta(m)*phi};
end for;

V := VectorSpace(F,d);

//ranks_of_3 := [Rank(sub<V | [V ! Eltseq(phi) : phi in subs]  >) : subs in Subsets(phis,3)];

ranks_of_4 := [Rank(sub<V | [V ! Eltseq(phi) : phi in subs]  >) : subs in Subsets(phis,4)];

//rank2_of_3 := [s : s in Subsets({m : m in M},3) | Rank(sub<V | [V ! Eltseq(Delta(m)*phi) : m in s]>) eq 2];
//#rank2_of_3;

function count(A)
    S := {x : x in A};
    c := {};
    for x in S do
	c := c join {[x,#[y : y in A | x eq y]]};
    end for;
    return c;
end function;

//count(ranks_of_3);
count(ranks_of_4);	
		   




      


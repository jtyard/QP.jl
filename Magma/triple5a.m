
function VtoP(v)
    m := Matrix([[c(v[i])*v[j] : j in [1..d]] : i in [1..d]]);
    return m/Trace(m);
end function;

P0 := VtoP(Eltseq(Transpose(phi)[1]));

ZZ := Integers();
ZZd := Integers(d);

"building SIC-POVM...";

PP := [[X^a*Z^b*P0*Z^-b*X^-a : b in [0..d-1]] : a in [0..d-1]];

P := func<j| PP[1 + ZZ!ZZd!j[1]][1+ZZ!ZZd!j[2]]>;

Tr3 := func<i,j,k|Trace(P(i)*P(j)*P(k))>;

Tr2 := func<i,j|Trace(P(i)*P(j))>;

//Strips := {Tr3(m, m*S[2], m*S[3]) : m in M};
//"computing all triple phases...";
//alltrips := [];
//SM := {m : m in M};
//n := 0;
//for s in Subsets(SM,3) do
//    ss := [m : m in s];
//    a := Tr3(ss[1],ss[2],ss[3]);
//    alltrips := alltrips cat [a];
//    n := n + 1;
//    if 100*Floor(n/100) eq n then
//	print n;
//    end if;
//end for;

//"counting distinct triple phases...";

//"found ", #{a : a in alltrips}, " triples";

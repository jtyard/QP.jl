S10 := SymmetricGroup(10);

//a5 := {@ 1,2,3,4,5 @};
J5 := Subsets({1..10},5);
J2 := Subsets({1..10},2);
J1 := {{j} : j in [1..10]};
#J5;

function tovec(a)
    v := Vector([0,0,0,0,0,0,0,0,0,0]);
    for i in a do
        v[i] := 1;
    end for;
    return v;
end function;

V5 := {tovec(a) : a in J5};
V2 := {tovec(a) : a in J2};
V1 := {tovec(a) : a in J1};
Vall := {tovec(a) : a in Subsets({1..10})};     

// T := Transpose(Matrix([Random(Vall),Random(Vall),Random(Vall)])); #{v*T : v in V1}; 

//TT := Subsets(V5,3);
/*
for t in T do
    if #{v*t : v in V5} eq 1024 then
        print t;
        break;
    end if;
end for;        
*/
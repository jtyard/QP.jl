S9 :=SymmetricGroup(9);
alpha := S9!(1,4,7)(2,5,8)(3,6,9);
gamma := S9!(1,2,3)(4,5,6)(7,8,9);

// beta := S9!(1,2,3)(4,6,5);
beta := S9!(4,5,6)(7,9,8);

minusone := S9!(1,9)(2,8)(3,7)(4,6);
tau := S9!(1,2)(4,5)(7,8);

G27 := sub<S9 | alpha,beta,gamma>;
G9  := sub<S9 | gamma, beta>;
H9  := sub<S9 | alpha,gamma>;
H18 := sub<S9 | alpha,gamma,minusone>;
G56 := sub<S9 | H18,beta>;
 
s221 := {1,2,4,5,9};
s32 := {1,2,4,5,6}; // fixed by beta, so really an H9 orbit is enough

simps := s221^G27 join s32^G27;

KC := SimplicialComplex(simps);

function star(p)
    return SimplicialComplex({s : s in Faces(KC,k), k in [1..5] | p in s});
end function;

function link(p)
    starp := star(p);
    return SimplicialComplex({s : s in Faces(starp,4) | not p in s});
end function;

function counts(C)
    return [#Faces(C,k) : k in [1..Dimension(C)+1]];
end function;


al := -1; load "3a.m";


//vtoj := [ [1 ,1], [2, 1], [0,1], [1,2], [2,2], [0,2], [1,0], [2,0], [0,0] ];

//jtov := [[9, 3, 6], [7, 1, 4], [8, 2, 5]];

jtov := [[1,2,3],[4,5,6],[7,8,9]];
vtoj := [ [0,0], [0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2]];

function sl2perm(s)
    return S9 ! [Index(vtoj,Eltseq((M ! vtoj[a])*Transpose(s))) : a in [1..9]];
end function;

function Tr9(s)
    ss := [a : a in s];
    return Tr3(vtoj[ss[1]],vtoj[ss[2]],vtoj[ss[3]]);
end function;

RTr9 := func<s|Tr9(s) + c(Tr9(s))>;

KC3 := Faces(KC,3);
KC4 := Faces(KC,4);
KC5 := Faces(KC,5);
KC6 := Faces(KC,6);

Rphases := {@ RTr9(t) : t in KC3 @};

phasecorbits := {@ {@ a, c(a) @} : a in phaseset @};


[#[s : s in KC4 | t subset s ] : t in KC3 | 8*Tr9(t) eq 1];
[#[s : s in KC4 | t subset s ] : t in KC3 | 8*Tr9(t) in {z3,z3^-1}];
[#[s : s in KC4 | t subset s ] : t in KC3 | 8*Tr9(t) in {z6,z6^-1}];
[#[s : s in KC4 | t subset s ] : t in KC3 | 8*Tr9(t) eq -1];

t5 := {@ t : t in KC3 | #[s : s in KC5 | t subset s ] eq 5 @};
t4 := {@ t : t in KC3 | #[s : s in KC5 | t subset s ] eq 4 @};
t3 := {@ t : t in KC3 | #[s : s in KC5 | t subset s ] eq 3 @};
t6 := {@ t : t in KC3 | #[s : s in KC5 | t subset s ] eq 6 @};


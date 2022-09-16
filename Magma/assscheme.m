// Given a set of vectors, do their inner products determine an association scheme?

//VV := [z3^k * Delta(j)*phi : j in M, k in [0..2]];
VV := [(z6)^k*Delta(j)*phi : j in M, k in [0]];

ips := {@ ip(v,w)/phinorm: v, w in VV @};

Ra := [[[v,w] : v,w in VV | ip(v,w)/phinorm eq a] : a in ips];

[#R : R in Ra];


for aa,bb,cc in [2..#ips] do
    counts := {@ #[z : z in VV | ip(xy[1],z)/phinorm eq ips[aa] and ip(z,xy[2])/phinorm eq ips[bb]] : xy in Ra[cc] @};
    if counts ne {@ 0 @} then
	print aa,bb,cc ,#Ra[aa], #Ra[bb], #Ra[cc],  counts;
    end if;
end for;

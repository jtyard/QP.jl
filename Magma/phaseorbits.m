// requires 5a.m

// "Computing triple phases";
//phases := [[[[ Tr3([0,0], [j1,j2], [k1,k2]) : k2 in [0..d-1]] : k1 in [0..d-1]] : j2 in [0..d-1]] : j1 in [0..d-1]];

//phaseset := {@ phases[j1,j2,k1,k2] : j1,j2,k1,k2 in [1..d] | [j1,j2] ne [1,1] and [k1,k2] ne [1,1] and [j1,j2] ne [k1,k2] @};

pstemp := phaseset;
sigsc := [c^j*sig^k : k in [0..7], j in [0,1]];
orbits := [];
while pstemp ne {} do
    a := pstemp[1];
    orbit := {@ sigsc[k](a) : k in [1..16] @};
    Append(~orbits,orbit);
    pstemp := pstemp diff orbit;
end while;

sigs := [sig^k : k in [0..7]];
realset := {@ (a+c(a))/2 : a in phaseset @};
rstemp := realset;
rorbits := [];
while rstemp ne {} do
    r:= rstemp[1];
    orbit := {@ sigs[k](r) : k in [1..8] @};
    Append(~rorbits,orbit);
    rstemp := rstemp diff orbit;
end while;

function phasetopairs(a)
    return [[M ! [j1,j2],M! [k1,k2]] : j1, j2, k1, k2 in [0..d-1] | phases[j1+1][j2+1][k1+1][k2+1] eq a];
end function;


    

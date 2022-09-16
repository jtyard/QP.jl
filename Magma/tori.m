gl2 := GL(2,Integers(d));
sl2 := SL(2,Integers(d));

alltori := AbelianSubgroups(gl2);

index :=  [i : i in [1..#alltori] | alltori[i]`order eq d^2 - 1];
torus := alltori[index[1]];
nsTgl := [T : T in Class(gl2,torus`subgroup)];
nsTsl := [T meet sl2 : T in nsTgl];

index :=  [i : i in [1..#alltori] | alltori[i]`order eq (d-1)^2];
torus := alltori[index[1]];
sTgl := [T : T in Class(gl2,torus`subgroup)];
sTsl := [T meet sl2 : T in sTgl];

#nsTgl, #sTgl;


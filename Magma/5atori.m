tori := AbelianSubgroups(gl2);

maxsplit_torus := tori[18]`subgroup;

maxnsplit_torus := tori[19]`subgroup;
		
diag := {gl2 ! [1, 0, 0, a] : a in [1,2,3,4]};

splitdiag := [T : T in Class(gl2,maxsplit_torus) | diag subset {t : t in T}];
nsplitdiag := [T : T in Class(gl2,maxnsplit_torus) | diag subset {t : t in T}];

#Class(gl2,maxsplit_torus);
#Class(gl2,maxnsplit_torus);
#splitdiag;
#nsplitdiag;


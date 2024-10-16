## Class fields
https://www.thofma.com/Hecke.jl/stable/

Some relevant functions from Hecke:
- Relative automorphism generators from  `Hecke.automorphism_groupQQ`
- Defining ray class fields: `rcf = ray_class_field((iseven(d) ? 2*d : d)*OK,infinite_places(K))` for `K = quadratic_field((d-3)*(d+1))[1]` and `OK = maximal_order(K)`.
- `number_field(rcf)` computes generators `absolute_automorphism_group(rcf)` of the full automorphism group and now works for all `d`.   
- TO DO: check that we can construct the ring ray class fields.
- `MapClassGrp` : {quotient of the class group} -> {ideals} 
- `MapRayClassGrp` : {quotient of a ray class group} -> {ideals prime to the conductor}
- `ClassField` 
- `ClassField_pp` Cyclic class field of prime-power degree
- `artin_map(rcf)` gives map from the ideal group to the set of automorphisms of `number_field(rcf)` i.e. from a `FacElemMon{Hecke.NfAbsOrdIdlSet{AnticNumberField, nf_elem}}` to a `Hecke.NfMorSet{NfRelNS{nf_elem}}`
- `complex_conjugation(F,infplace)` extends the Artin map to infinite places, giving complex conjugation in the corresponding complex embeddings. 
- `automorphism_group(rcf)` gives a map from a `GrpGen` to the set of automorphisms of `number_field(rcf)` fixing the base, and `inv(rcf.quotientmap)` works. 
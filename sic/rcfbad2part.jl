using Oscar 
using QP

#d=7

# Okay so absolute_automorphism_group fails if the field is not Galois over QQ - that's the issue.
# But what other problems do we have 




set_verbose_level(:ClassField, 1)

function ab_aut(d)
    K = quadratic_field((d-3)*(d+1))[1]
    OK = maximal_order(K)
    rcf = ray_class_field((isodd(d) ? d : 2*d)*OK, real_places(K))
    absolute_automorphism_group(rcf)
end

K = quadratic_field(3*5*2)[1]
OK = maximal_order(K)
rcf = ray_class_field(27*OK, real_places(K))
absolute_automorphism_group(rcf)


# Without infinite places, fails for d = 16,20 but d=27 runs forever.

#rcf1 = ray_class_field((isodd(d) ? d : 2*d)*OK,[real_places(K)[1]])
#rcf2 = ray_class_field((isodd(d) ? d : 2*d)*OK,[real_places(K)[2]])
#rcf0 = ray_class_field((isodd(d) ? d : 2*d)*OK)
#number_field(rcf,using_stark_units=true);

# absolute_automorphism_group(SICdata(d).rcf)

# All the even dimensions below computed the wrong thing (d vs)
# works for d in [4,6,7,8,10,12,14,18,19,35]

# Otherwise gets caught lifting automorphisms of degree 2^n cyclic class field:
# Computing rel_auto 0.0032291 @ aut_L_rel = rel_auto(C)::Vector{NfRelNSToNfRelNSMor_nf_elem}
# Extending automorphismsExtending auto pp ...

# runs forever for d in [5,9,11,13,15,17,21] (deg [16, 72, 80, 192, 48,192,192] 2-parts [2^4, 2^3, 2^4, 2^6, 2^4, 2^6, 2^6])
# For d in [16,20] we get (deg [256,384] 2parts [2^8, 2^7])
# ERROR: LoadError: Not yet implemented

# Weird.

# Oops I think this was with K = Q(sqrt(2))
# works for d in [4,6,7,12,14,19,35]
# runs forever for d in  [5,9,10,11,13,15,17,21]
# not implemented error for d in [8,10,16,18,20]
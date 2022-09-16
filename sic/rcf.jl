using Oscar 

set_verbose_level(:ClassField, 1)
K = quadratic_field(2)[1]
OK = maximal_order(K)
rcf = ray_class_field(7*OK,real_places(K))
#number_field(rcf,using_stark_units=true);
absolute_automorphism_group(rcf);


# absolute_automorphism_group(SICdata(d).rcf) 
# works for d in [4,6,7,12,19,35]
# runs forever for d in  [5,9,10,11,13,15,17,21]
# not implemented error for d in [8,10,14,16,18,20]
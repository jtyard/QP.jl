

# Given a number field F/K, returns a simple number field over K and an isomorphism to F/K
function relative_simple_field(F)
    K = base_field(F)
    FQ,ph = absolute_simple_field(F)
    L,psi = Hecke.subfield(FQ,[preimage(ph,F(gen(K)))])
    relative_simple_extension(FQ,L)
end



# N = 5
#S = SicData(N)
#rcf = ray_class_field((isodd(S.d) ? S.d : 2*S.d)*S.OK,S.inf)
#@time F = number_field(rcf,using_stark_units=true)

#@time is_isomorphic(absolute_simple_field(F)[1],absolute_simple_field(S.F)[1])


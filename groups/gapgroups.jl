# Figuring out how to build groups in Oscar

# GAP.Packages.install("forms"); (Didn't have to do this any more)
# ] add "GAP"
# GAP.LoadPackageAndExposeGlobals("forms","Forms")

@time using Oscar

gap = GAP.Globals

GAP.LoadPackageAndExposeGlobals("forms","Forms")

quadratic_form = Forms.QuadraticFormByMatrix

z0 = gap.Z(2)*0
z1 = gap.Z(2)^0

mat = GAP.julia_to_gap([z0 z1; z0 z0])

q = quadratic_form(mat)

#  GAP.Globals.OrbitStabilizer 
#  https://www.gap-system.org/Manuals/doc/ref/chap41.html#X797BD60E7ACEF1B1


# Hopefully they fix this - is there a more general type?
import Oscar.transpose
transpose(M::Main.ForeignGAP.MPtr) = GAP.Globals.TransposedMatImmutable(M)

# Wow I can make the Kuhnel CP^2 with 
# Polymake.topaz.complex_projective_plane()
# And RP^2 is 
# Polymake.topaz.real_projective_plane()
# ∑ |n⟩|n⟩/√n


#n = 10  # 10 x 10 matrices: TODO: test with a range of different degrees
#
#testdata = []
#
#for (p,d) in [(2,1), (17,1), (17,2)]
#  q = p^d
#  F, _ = FiniteField(ZZRingElem(p), d, "a")
#  mats = []
#  push!(mats, (one(GL(n,q)), "GAP GF("*string(p)*","*string(d)*")"))
#  push!(mats, (one(GL(n,q)).X, "raw GAP GF("*string(p)*","*string(d)*")"))
#  push!(mats, (identity_matrix(F, n), "FiniteField("*string(p)*","*string(d)*")"))
#  if d == 1
#    push!(mats, (identity_matrix(GF(p), n), "Oscar GF("*string(p)*","*string(d)*")"))
#  end
#  push!(testdata, (mats, "GF("*string(p)*","*string(d)*")"))
#end
export SicData

mutable struct SicData
    d::ZZRingElem
    D::ZZRingElem    
    D0::ZZRingElem
    f::ZZRingElem
    P::MatPolyRing
    K::NumField
    inf::Vector{InfPlc}
    OK::NumFieldOrder # the maximal order Z[uf]
    uf::NumFieldOrderElem #or NumFieldElem?
    OD::NumFieldOrder # the minimal order Z[b]
    Ob::NumFieldOrder # the minimal order Z[b]
    b::NumFieldOrderElem # or NumFieldOrdElem?
    rcf::ClassField
    F::NumField
    c::NumFieldHom
    e::Hecke.NumFieldEmb # or HeckeNumFieldEmbNfNS ?
    I::MPolyIdeal
    Iall::MPolyIdeal
    #ring_class_field::ClassField
    #function SicData(d::Int;build_nf=true)
    @memoize function SicData(d::Int;build_nf=true)
        if d == 2
            P = matrix_polynomial_ring(QQ,d)
            S = new(2,-3,-3,1,P,rationals_as_number_field()[1])
            S. I = Ihplus(P) + Iminors(P,2)
            S. Iall = Ihplus(P) + Ihminus(P)
            return S
        end
        if d == 3
            P = matrix_polynomial_ring(QQ,d)
            S = new(3,0,0,1,P,rationals_as_number_field()[1])
            S.I = Ihplus(P) + Iminors(P,2)
            S.Iall = Ihplus(P) + Ihminus(P)
            return S
        end
        D = ZZ((d-3)*(d+1))
        D0 = fundamental_discriminant(D) 
        f = ZZ(sqrt(D//D0)) 
        #println("making projective space")
        P = matrix_polynomial_ring(QQ,d)
        #println("Ihplus + Iminors")
        #I = Ihplus(P) + Iminors(P,2)
        #println("Ihplus + Ihminus")
        #Iall = Ihplus(P) + Ihminus(P)
        K = quadratic_field(Hecke.squarefree_part(D0))[1] 
        #println("finding real places")
        inf = real_places(K)
        #println("computing maximal order")
        OK = maximal_order(K)
        #println("finding fundamental unit")
        uf = d > 3 ? fundamental_unit(OK) : OK(1)
        OD = quadratic_order(D)
        bb = (d-1 + sqrt(K(D)))//K(2)
        Ob = quadratic_order(bb)
        b = Ob(bb)
        #println("creating rcf")
        rcf = ray_class_field((isodd(d) ? d : 2*d)*OK,infinite_places(K))
        #println("created rcf")
        #println("Constructing new")
        S = new(
        d,
        D,
        D0,
        f,
        P,
        K,
        inf,
        OK,
        uf,
        OD,
        Ob,
        b,
        rcf)
        if build_nf
            println("Constructing number field")
            S.F = number_field(rcf)
            #F = S.F
            #S.F = number_field(rcf,using_stark_units = true)
            #println("Finding complex conjugation")
            println("Finding complex conjugation")
            S.c = complex_conjugation(rcf,inf[2])
            
            # would be nice if we could do
            # complex_conjugation(S.F) = c 
            # but it's not so easily possible

            #CC = AcbField(64) #bitsize should be increased for larger fiducials but this will work for small ones. Revisit this.
            #zdd = exp(2*pi*onei(CC)/(2*d))
            #println("Computing embeddings")
            #embs = complex_embeddings(F)
            #println("Finding an embedding")
            # Next line scales badly and might not even be necessary
            #S.e = [e for e in embs if overlaps(e(zetaN(2*d,F)),zdd) && overlaps(e(F(gen(S.K))),sqrt(CC(Hecke.squarefree_part(D0))))][1]
        end
        #println("Returning structure")
        S
    end
end
#First julia code I ever wrote, wraps pari to get special values of L-functions relevant for Stark units.  Much more stable pure julia options already in Oscar
using PyCall

# 7*, 13, 19*, 31, 37, 43, 61, 67*, 73, 79, 97, 103*, 109, 127, 139, 151, 157, 163, 181, 193, 199*,...
# * means hplus = h. All d* in this list have h=1, except for 103* which has h = 2.

# 
function lvals(d::Integer; verbose=false,rchi=false)
    #setprecision(64)
    cypari = pyimport("cypari")
    pari = cypari.pari
    PREC = precision(BigFloat)
    pari.set_real_precision_bits(PREC)   
    pari.allocatemem(10^9)

    D = (d-3)*(d+1)
    x = pari("x")
    K = pari.bnfinit(x^2-D, 1,precision=PREC)
    
    P = K.idealprimedec(d)
    P1 = get(P,0)
    P2 = get(P,1)

    H = K.bnrinit([1,[0,0]],1)
    Hinf = K.bnrinit([1,[1,1]],1)
    E1 = K.bnrinit([P2,[1,0]],1)
    E2 = K.bnrinit([P2,[0,1]],1)
    E12 = K.bnrinit([P2,[1,1]],1)
    
    n1 = E1.getattr("clgp").getattr("no").python()
    n2 = E2.getattr("clgp").getattr("no").python()
    
    if n1 > n2
        P1,P2 = P2,P1
        E1,E2 = E2,E1
        n1,n2 = n2,n1
    end

    C = H.getattr("cyc")
    Cinf = Hinf.getattr("cyc")
    println("## d = ", d,C==Cinf ? "*" : "")
    #println(E12.getattr("cyc"))
    #println(E1.getattr("cyc"))
    println(E2.getattr("cyc")," ", n2," ",n1)   
    println(Cinf) 
    println(C)
    println()

    bnrL1 = E2.bnrL1(0,6,precision=PREC)
    
    Lpy = [get(get(get(bnrL1,i),1),1) for i = 0:bnrL1.length()-1]
    L = [convert(BigFloat,a.real()) + convert(BigFloat,a.imag())*im for a in Lpy]

    R = [convert(Integer,get(get(get(bnrL1,i),1),0)) for i = 0:bnrL1.length()-1]
    Chi = [convert(Array{Integer},get(get(bnrL1,i),0)) for i = 0:bnrL1.length()-1]

    if verbose
        println(E2.getattr("nf"))
        println()
        println(L)
        println()
        println(R)
        println()
        println(Chi)
        println()
    end

    if rchi
        L,R,Chi
    else 
        L
    end
end
# PyObject [[[1], [1, 0.632974319]], [[0], [2, -0.857536904]]]
# vs 
# PyObject [[1, 0.632974319], [2, -0.857536904]]

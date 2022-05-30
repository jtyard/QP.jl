
using PyCall
cypari = pyimport("cypari")
pari = cypari.pari
# 7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109, 127, 139, 151, 157, 163, 181, 193, 199

#function lvals(d::Integer)
for d in [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109, 127, 139, 151, 157, 163,181, 193, 199]
    D = (d-3)*(d+1)
    x = pari("x")
    K = pari.bnfinit(x^2-D)
    
    P = K.idealprimedec(d)
    P1 = get(P,0)
    P2 = get(P,1)

    H = K.bnrinit([1,[0,0]],1)
    Hinf = K.bnrinit([1,[1,1]],1)
    E1 = K.bnrinit([P2,[1,0]],1)
    E2 = K.bnrinit([P2,[0,1]],1)
    E12 = K.bnrinit([P2,[1,1]],1)
    
    C = H.getattr("cyc")
    Cinf = Hinf.getattr("cyc")

    println("## d = ", d,C==Cinf ? "*" : "")
    println(E12.getattr("cyc"))
    println(E1.getattr("cyc"))
    println(E2.getattr("cyc"))   
    println(Cinf) 
    println(C)
    println()
end

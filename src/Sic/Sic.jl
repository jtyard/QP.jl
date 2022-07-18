using Oscar
using Caching


# These are all terrible names to be global in QP.jl
# ultimately will move them out of global after things are working
export minors, Im, QQXp
export XX, h, h2, hp, hm
export Ih, Ih2, Ihm, Ihp, Ic, Icc, Itr0, IT


# Making new rings for the overlaps - maybe use routines from Vars.jl instead? 
#Aij(i::Union{Int,nmod},j::Union{Int,nmod},N::Int) = gens(CA(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
#Aij(j::nmod_mat) = Aij(j[1],j[2],Int(characteristic(base_ring(j))))
#Aij(N::Int) = matrix(CA(N),[[Aij(i,j,N) for j =0:N-1] for i =0:N-1])
#
#aij(j::nmod_mat) = tr(Aij(Int(characteristic(base_ring(j))))*heis(j))

X = Xij

function XX(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    sum([X(a,a+j[1],N)*X(a+j[1]+j[2],a+j[2],N) for a in 0:N-1])
end

XX(j1,j2,N::Int) = XX(ZN(N)[j1 j2])

function h(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    (N+1)*XX(j) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*TrX(N)^2
end
h(j1,j2,N::Int) = h(ZN(N)[j1 j2])

function h2(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    (N+1)*XX(j) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*TrX2(N)
end
h2(j1,j2,N::Int) = h2(ZN(N)[j1 j2])

function hp(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    ((N+1)//2)*(XX(j) + XX(ZN(N)[j[2] j[1]])) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*(TrX(N)^2 + TrX2(N))*(1//2)
end
hp(j1,j2,N::Int) = hp(ZN(N)[j1 j2])

function hm(j::nmod_mat) 
    N = Int(characteristic(base_ring(j)))
    ((N+1)//2)*(XX(j) - XX(ZN(N)[j[2] j[1]])) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*(TrX(N)^2 - TrX2(N))*(1//2)
end
hm(j1,j2,N::Int) = hm(ZN(N)[j1 j2])

############
# Useful ideals
############

Ih(N::Int) = ideal(QQX(N),[h(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ih2(N::Int) = ideal(QQX(N),[h2(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihp(N::Int) = ideal(QQX(N),[hp(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihm(N::Int) = ideal(QQX(N),[hm(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])

Im(N::Int) = ideal(QQX(N),minors(N))

Ic(N::Int) = ideal(QQX(N), [Xij(j1,j2,N) - Xij(j2,j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
Icc(N::Int) = ideal(QQX(N), [Xij(j1,j2,N) - Xij(-j2,-j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
#Itr1(N) = ideal([TrX(N) - 1])
Itr0(N::Int) = ideal(QQX(N),[TrX(N)])
IT(a,N::Int) = ideal(QQX(N), [Xij(j1,j2) - Xij(a*j1,a*j2) for j1 in 0:N-1 for j2 in 0:N-1] )

# duh Itr0(N) = R(X_{0,0} + X_{1,1} + ⋯ ) subset RXp(N) = R({X_{i,j}})





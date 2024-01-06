
# Making new rings for the overlaps - maybe use routines from Vars.jl instead? 
#Aij(i::Union{Int,zzModRingElem},j::Union{Int,zzModRingElem},N::Int) = gens(CA(N))[1 + (Int(j) % N) + N*(Int(i) % N)]
#Aij(j::zzModMatrix) = Aij(j[1],j[2],Int(characteristic(base_ring(j))))
#Aij(N::Int) = matrix(CA(N),[[Aij(i,j,N) for j =0:N-1] for i =0:N-1])
#
#aij(j::zzModMatrix) = tr(Aij(Int(characteristic(base_ring(j))))*heis(j))

# Not actually used anywhere, subsumed by SICdata



X = Xij

function XX(j::zzModMatrix) 
    N = Int(characteristic(base_ring(j)))
    sum([X(a,a+j[1],N)*X(a+j[1]+j[2],a+j[2],N) for a in 0:N-1])
end

XX(j1,j2,N::Int) = XX(ZN(N)[j1 j2])

function h(j::zzModMatrix) 
    N = Int(characteristic(base_ring(j)))
    (N+1)*XX(j) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*TrX(N)^2
end
h(j1,j2,N::Int) = h(ZN(N)[j1 j2])

function h2(j::zzModMatrix) 
    N = Int(characteristic(base_ring(j)))
    (N+1)*XX(j) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*TrX2(N)
end
h2(j1,j2,N::Int) = h2(ZN(N)[j1 j2])

function hp(j::zzModMatrix) 
    N = Int(characteristic(base_ring(j)))
    ((N+1)//2)*(XX(j) + XX(ZN(N)[j[2] j[1]])) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*(TrX(N)^2 + TrX2(N))*(1//2)
end
hp(j1,j2,N::Int) = hp(ZN(N)[j1 j2])

function hm(j::zzModMatrix) 
    N = Int(characteristic(base_ring(j)))
    ((N+1)//2)*(XX(j) - XX(ZN(N)[j[2] j[1]])) - ((j[1]==0 ? 1 : 0) + (j[2]==0 ? 1 : 0))*(TrX(N)^2 - TrX2(N))*(1//2)
end
hm(j1,j2,N::Int) = hm(ZN(N)[j1 j2])


############
# Useful ideals
############

Ih(N::Int; graded = false) = ideal(QQX(N, graded = graded),[h(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ih2(N::Int; graded = false) = ideal(QQX(N, graded = graded),[h2(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihp(N::Int; graded = false) = ideal(QQX(N, graded = graded),[hp(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])
Ihm(N::Int; graded = false) = ideal(QQX(N, graded = graded),[hm(j1,j2,N) for j1 in 0:N-1 for j2 in j1:N-1])

Im(N::Int; graded = false) = ideal(QQX(N, graded = graded),minors(N))

Ic(N::Int; graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(j2,j1,N) for j1 in 0:N-1 for j2 in 0:N-1 if j1 < j2] )
Icc(N::Int; graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(-j2,-j1,N) for j1 in 0:N-1 for j2 in 0:N-1] )
Itr1(N) = ideal([TrX(N) - 1])
Itr0(N::Int; graded = false) = ideal(QQX(N, graded = graded),[TrX(N)])
IT(a,N::Int; graded = false) = ideal(QQX(N, graded = graded), [Xij(j1,j2,N) - Xij(ZN(N)(a)*j1,ZN(N)(a)*j2,N) for j1 in 0:N-1 for j2 in 0:N-1] )

# duh Itr0(N) = R(X_{0,0} + X_{1,1} + ⋯ ) subset RXp(N) = R({X_{i,j}})

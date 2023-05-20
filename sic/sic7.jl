


# 
N=7; Zd = ZN(N); s = Zd[0 -1; 1 0]; t = Zd[1 1; 0 1]; X = gpX(N); Z = gpZ(N);

#j = Zd[rand(0:N-1); rand(0:N-1)]; k = Zd[rand(0:N-1); rand(0:N-1)]; heis(j)*heis(k) ==  heis(j+k)*heiscocycle(j,k)
#heis(j)*heis(k)*heis(-j)heis(-k) == heispairing(j,k)*weil_U(Zd[1 0;0 1])

#j = Zd[rand(0:N-1); rand(0:N-1)]; weil_U(s)*heis(j)*weil_U(s)^-1 == heis(s*j)
#j = Zd[rand(0:N-1); rand(0:N-1)]; weil_U(t)*heis(j)*weil_U(t)^-1 == heis(t*j)


S = SicData(7)
K = S.K
OK = S.OK

s = sqrt(OK(2))

I = ideal(OK,2+s)
#I = ideal(OK,7)
rcf0 = ray_class_field(I)
rcf1 = ray_class_field(I,[real_places(K)[1]])
rcf2 = ray_class_field(I,[real_places(K)[2]])
rcf = ray_class_field(I,real_places(K))



#F = number_field(rcf)
#A = Hecke.rel_auto(rcf)


#@time Isat = saturation(I,Rp)
#hilbert_series_reduced(quo(R,I)[1])

#@time pd = primary_decomposition(I,alg = :SY)
#pp = [p[2] for p in pd]
#display([[dim(p+q) for q in pp] for p in pp])



# Calculations verify the following: 

# Irad = radical(I)
# Isat = saturation(I,Rp)
# pd = primary_decomposition(I)

# N=2 I = Irad = Isat; length(pd) = 2; dim = 0

# N=3 I = Isat neq Irad; length(pd) = 8; dim = 1
# julia> hilbert_series_reduced(quo(R,I)[1])
# (t^5 + 7*t^4 + 16*t^3 + 16*t^2 + 7*t + 1, t^2 - 2*t + 1)
# julia> hilbert_series_reduced(quo(R,radical(I))[1])
#  0.279213 seconds (169.81 k allocations: 6.843 MiB, 7.30% compilation time)
# (16*t^2 + 7*t + 1, t^2 - 2*t + 1)

# N=4 I = Isat; length(pd) = ?? (should be 56); dim = 0

#julia> @time Isat = saturation(I,Rp)
#754.226843 seconds (5.61 G allocations: 103.652 GiB, 0.19% gc time, 0.00% compilation time)
#ideal(X_{0,0}*X_{1,1} - X_{0,1}*X_{1,0}, X_{0,0}*X_{1,2} - X_{0,2}*X_{1,0}, X_{0,0}*X_{1,3} - X_{0,3}*X_{1,0}, X_{0,1}*X_{1,2} - X_{0,2}*X_{1,1}, X_{0,1}*X_{1,3} - X_{0,3}*X_{1,1}, X_{0,2}*X_{1,3} - X_{0,3}*X_{1,2}, X_{0,0}*X_{2,1} - X_{0,1}*X_{2,0}, X_{0,0}*X_{2,2} - X_{0,2}*X_{2,0}, X_{0,0}*X_{2,3} - X_{0,3}*X_{2,0}, X_{0,1}*X_{2,2} - X_{0,2}*X_{2,1}, X_{0,1}*X_{2,3} - X_{0,3}*X_{2,1}, X_{0,2}*X_{2,3} - X_{0,3}*X_{2,2}, X_{0,0}*X_{3,1} - X_{0,1}*X_{3,0}, X_{0,0}*X_{3,2} - X_{0,2}*X_{3,0}, X_{0,0}*X_{3,3} - X_{0,3}*X_{3,0}, X_{0,1}*X_{3,2} - X_{0,2}*X_{3,1}, X_{0,1}*X_{3,3} - X_{0,3}*X_{3,1}, X_{0,2}*X_{3,3} - X_{0,3}*X_{3,2}, X_{1,0}*X_{2,1} - X_{1,1}*X_{2,0}, X_{1,0}*X_{2,2} - X_{1,2}*X_{2,0}, X_{1,0}*X_{2,3} - X_{1,3}*X_{2,0}, X_{1,1}*X_{2,2} - X_{1,2}*X_{2,1}, X_{1,1}*X_{2,3} - X_{1,3}*X_{2,1}, X_{1,2}*X_{2,3} - X_{1,3}*X_{2,2}, X_{1,0}*X_{3,1} - X_{1,1}*X_{3,0}, X_{1,0}*X_{3,2} - X_{1,2}*X_{3,0}, X_{1,0}*X_{3,3} - X_{1,3}*X_{3,0}, X_{1,1}*X_{3,2} - X_{1,2}*X_{3,1}, X_{1,1}*X_{3,3} - X_{1,3}*X_{3,1}, X_{1,2}*X_{3,3} - X_{1,3}*X_{3,2}, X_{2,0}*X_{3,1} - X_{2,1}*X_{3,0}, X_{2,0}*X_{3,2} - X_{2,2}*X_{3,0}, X_{2,0}*X_{3,3} - X_{2,3}*X_{3,0}, X_{2,1}*X_{3,2} - X_{2,2}*X_{3,1}, X_{2,1}*X_{3,3} - X_{2,3}*X_{3,1}, X_{2,2}*X_{3,3} - X_{2,3}*X_{3,2}, 3*X_{0,0}^2 - 4*X_{0,0}*X_{1,1} - 4*X_{0,0}*X_{2,2} - 4*X_{0,0}*X_{3,3} + 3*X_{1,1}^2 - 4*X_{1,1}*X_{2,2} - 4*X_{1,1}*X_{3,3} + 3*X_{2,2}^2 - 4*X_{2,2}*X_{3,3} + 3*X_{3,3}^2, -X_{0,0}^2 + 3*X_{0,0}*X_{1,1} - 2*X_{0,0}*X_{2,2} + 3*X_{0,0}*X_{3,3} - X_{1,1}^2 + 3*X_{1,1}*X_{2,2} - 2*X_{1,1}*X_{3,3} - X_{2,2}^2 + 3*X_{2,2}*X_{3,3} - X_{3,3}^2, -X_{0,0}^2 - 2*X_{0,0}*X_{1,1} + 8*X_{0,0}*X_{2,2} - 2*X_{0,0}*X_{3,3} - X_{1,1}^2 - 2*X_{1,1}*X_{2,2} + 8*X_{1,1}*X_{3,3} - X_{2,2}^2 - 2*X_{2,2}*X_{3,3} - X_{3,3}^2, 5*X_{0,1}*X_{2,1} + 5*X_{0,3}*X_{2,3} + 5*X_{1,0}*X_{3,0} + 5*X_{1,2}*X_{3,2}, 5*X_{0,1}*X_{3,2} + 5*X_{0,3}*X_{1,2} + 5*X_{1,0}*X_{2,3} + 5*X_{2,1}*X_{3,0}, 5*X_{0,1}*X_{0,3} + 5*X_{1,0}*X_{1,2} + 5*X_{2,1}*X_{2,3} + 5*X_{3,0}*X_{3,2}, 5*X_{0,2}^2 + 5*X_{1,3}^2 + 5*X_{2,0}^2 + 5*X_{3,1}^2, 5*X_{0,2}*X_{1,3} + 5*X_{0,2}*X_{3,1} + 5*X_{1,3}*X_{2,0} + 5*X_{2,0}*X_{3,1})

#julia> I == Isat
#true



#dim(affine_cone(S)) - 1


#@time pd = primary_decomposition(I,alg = :SY)
#pp = [p[2] for p in pd]
#display([[dim(p+q) for q in pp] for p in pp])



#Isat = saturation(I,Rp)
#pdsat = primary_decomposition(Isat)
#A,_ = quo(RX(N),Isat + Itr)

#G1 = abelian_group([0])
#R1 = grade(RX(N))

#G12 = abelian_group([0,N])
#R12 = grade(RX(N),[G12([1,i-j]) for i in 0:N-1 for j in 0:N-1])


#B,_ = quo(RG,ideal([RG(e) for e in eqns]))

#dims = [[dim(homogenous_component(R12,G12([a,b]))[1]) for b in 0:N-1] for a in 0:4]

#println(Array(dims))
#println([sum(d) for d in dims])

#println([dim(homogenous_component(R1,G1([a]))[1]) for a in 0:3])

#S,Yvars = PolynomialRing(QQ,[string("Y_{",i,",",j,"}") for i in 0:N-1 for j in 0:N-1])
#Y(j) = Yvars[1 + mod(j[2],N) + N*mod(j[1],N)]
#Y(j1,j2) = Y([j1,j2])
#G123 = abelian_group([0,N,N])
#S123 = grade(S,[G123([1,i,j]) for i in 0:N-1 for j in 0:N-1])
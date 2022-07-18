using Oscar
using Plots
#unicodeplots()
gr()
#pgfplotsx()

size=10
param = (1.0/size):.01:size
X = -size:.1:size

# Axes
plot(xlims = (-size,size), xticks = -size:size, 
    ylims = (-size,size), yticks = -size:size, 
    aspect_ratio=:equal, 
    draw_arrow=:arrow, 
    framestyle=:origin,
    legend=false
    )
# hyperbolas
plot!(param,x->x^-1, color=:black)
plot!(-param,x->x^-1, color=:black)
plot!(param,x->-x^-1, color=:red)
plot!(-param,x->-x^-1, color=:red)

# diagonal
plot!(X,x-> x, color=:black)
#plot!(X,x->-x)

d = 7
K,s2 = quadratic_field(2)
OK = maximal_order(K)

# uf    =  1 + s2 
# uf^-1 = -1 + s2
# up = uf^2


embed(x) = [v(x) for v in real_places(parent(x))]

uf = 1+sqrt(2)

b = 1/uf^2
binv = uf^2

O_D0 = [m*Vector([1,1]) + n*Vector([-1/uf,uf]) for m in -8:8 for n in -8:8]
sublattice1 = [m*Vector([uf,-1/uf]) + n*Vector([-1/uf,uf]) for m in -4:4 for n in -4:4] # sqrt(2) OK
sublattice2 = [m*Vector([uf^2,1/uf^2]) + n*Vector([1/uf^2,uf^2]) for m in -3:3 for n in -3:3] # sqrt(2)
#D = (sqrt(32),-sqrt(32))
plot!(seriestype=:scatter,[Tuple(p) for p in O_D0],color=:white)
plot!(seriestype=:scatter,[Tuple(p) for p in sublattice1],color=:orange)
plot!(seriestype=:scatter,[Tuple(p) for p in sublattice2],color=:black)



#plot!(seriestype=:scatter,[c1,c2],color=:blue)

#plot!(seriestype=:scatter,[(uf,-1/uf),(1/uf,-uf),(-1/uf,uf),(-uf,1/uf)], color=:red)
#plot!(seriestype=:scatter,[(uf^2,1/uf^2),(1/uf^2,uf^2),(-uf^2,-1/uf^2),(-1/uf^2,-uf^2)], color=:black)
#plot!(seriestype=:scatter,[(0,0),(6,6),(-6,-6)],color=:black)
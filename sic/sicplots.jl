using Plots
#unicodeplots()
gr()
#pgfplotsx()

units = (1.0/7):.01:7
X = -7:.1:7

plot(xlims = (-7,7), xticks = -7:7, 
    ylims = (-7,7), yticks = -7:7, 
    aspect_ratio=:equal, 
    draw_arrow=:arrow, 
    framestyle=:origin,
    legend=false
    )
plot!(units,x->x^-1, color=:blue)
plot!(-units,x->x^-1, color=:blue)
plot!(units,x->-x^-1, color=:red)
plot!(-units,x->-x^-1, color=:red)
plot!(X,x-> x, color=:black)
#plot!(X,x->-x)

plot!(seriestype=:scatter,[(1,1),(-1,-1)],color=:purple)

c1 = (sqrt(3),-sqrt(3))
c2 = (3+2*sqrt(3), 3-2*sqrt(3))
plot!([c1,c2],color=:blue)
plot!(seriestype=:scatter,[c1,c2],color=:blue)

u = 2+sqrt(3)
plot!(seriestype=:scatter,[(1,1),(-1,-1),(u,1/u), (1/u,u),(-u,-1/u), (-1/u,-u)],color=:blue)

uf = 1+sqrt(2)
plot!(seriestype=:scatter,[(uf,-1/uf),(1/uf,-uf),(-1/uf,uf),(-uf,1/uf),
(uf^2,1/uf^2),(1/uf^2,uf^2),(-uf^2,-1/uf^2),(-1/uf^2,-uf^2)], color=:red)
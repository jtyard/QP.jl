# zonal polynomial for Johnson scheme
# following SPLAG p. 256 (22)

johnson_zonal(n,k,t) = sum([(-1)^i * binomial(k,i)*binomial(n+1-k,i)*binomial(t,i) // (binomial(w,i)*binomial(n-w,i))])


#J(n,k) = Polymake.graph.johnson_graph(n,k)
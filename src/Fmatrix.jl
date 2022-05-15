using Oscar

_, q = RationalFunctionField(QQ,"q")

qint(n) = sum([q^(n-1 - 2*i) for i in 0:n-1])

zd(d) = cyclotomic_field(d)[2]

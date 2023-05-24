using Oscar

G = Sp(4,2)

f1 = G.ring(1)
f0 = G.ring(0)

B = kron([f1 f0; f0 f1],[f0 f1; f1 f0])
B2 = kron([f0 f1; f1 f0],[f1 f0; f0 f1])

#gap = GAP.Globals

#GAP.LoadPackageAndExposeGlobals("forms","Forms")

# maxhorn on Slack says to use 
# GAP.Packages.load("forms"; install=true)
# which returns True though I wonder if it could already have been loaded

#qf = Forms.QuadraticFormByMatrix
#ζ×

# Probably already in some Sage notebook
# Known (relative) class numbers for 2^nth cyclotomic with n = 3:9
hm = [1,1,1, 17, 359057, 10449592865393414737, 6262503984490932358745721482528922841978219389975605329] 
# For comparison: the main bound
hmm = [BigInt(2)^(n*2^(n-3)) for n in 3:9]
export hm, hmm

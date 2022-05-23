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
#ζ
using Symbolics, SymbolicUtils
using SymbolicUtils.Rewriters

using Symbolics, SymbolicUtils

@variables t, A(t), x, y, z

# 1. Create a raw differential term
D = Differential(t)
D2 = Differential(t)^2 
# This creates a structure: Differential(t)(A(t))
target = D(A) + 5*A

# 2. Define a rule to replace D(A) with 0
rule = @acrule((D(D(~a))/D(D(~b))) => 0.0)
rule2 = Postwalk(rule)

r

simplify((D2(A)/D2(x)*(4x-z)+ y), rewriter = rule2, expand=true)

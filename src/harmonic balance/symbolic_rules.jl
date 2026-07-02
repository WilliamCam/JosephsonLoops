using SymbolicUtils.Rewriters

# 2. Define a rule to replace D(t) with 0
remove_slow_d2_rule_1 = @acrule (~c * D(~a)) => 0.0
remove_slow_d2_rule_2 = @acrule (D(~a)) => 0.0
remove_slow_d2_rewriter = [remove_slow_d2_rule_1, remove_slow_d2_rule_2] |> Chain |> Postwalk

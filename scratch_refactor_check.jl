# Verify the stored-`vars` refactor: linearised_vars returns the jacobian column order,
# linearised_row_map derives from it (reverse-lookup in variable_map), and get_output still
# reads the right rows of the linearised response.
using JosephsonLoops
const jls = JosephsonLoops
using Symbolics, ModelingToolkit, Printf

loops = [["P1", "C1", "J1"]]
circuit = jls.process_netlist(loops)
model, u0, guesses = jls.build_circuit(circuit)
sys = jls.HarmonicSystem(model, jls.P1.Isrc.ω, 2; determine_jacobian=true)

# 1. jacobian_vars is stored and sized to the jacobian
lv = jls.linearised_vars(sys)
@printf("linearised_vars length = %d, jacobian dim = %d, match = %s\n",
        length(lv), size(sys.jacobian[1], 1), length(lv) == size(sys.jacobian[1], 1))

# 2. row map derived by reverse-lookup — print it, and check it's a bijection onto 1:n
rm = jls.linearised_row_map(sys)
println("\nrow map (sorted by index):")
for k in sort(collect(keys(rm)), by = x -> rm[x]); println(lpad(rm[k], 3), "  ", k); end
@printf("\nrow map covers every column once: %s\n",
        sort(collect(values(rm))) == collect(1:length(lv)))

# 3. consistency: column idx's symbol == variable_map[key] that maps to idx
consistent = all(isequal(Symbolics.unwrap(lv[idx]), Symbolics.unwrap(sys.variable_map[key]))
                 for (key, idx) in rm)
@printf("vars[idx] == variable_map[key] for all rows: %s\n", consistent)

# 4. end-to-end: get_output must equal a manual row read of the response
ps = Dict(jls.P1.Isrc.I => 0.0, jls.P1.Rsrc.R => 50.0, jls.C1.C => 100.0e-15,
    jls.J1.C => 1000.0e-15, jls.J1.I0 => jls.Φ₀/(2π*1000.0e-12), jls.J1.R => 10e3)
ω_vec = collect(2*pi*(4.75:0.005:4.85)*1e9)
U = zeros(Float64, length(lv)); U[rm[("P1₊i", 1, :Sin)]] = 1.0e-9     # manual unit kick
lin = jls.HarmonicProblem(sys, ω_vec, ps; linear_response=(2*pi*4.80e9, U))
jls.solve!(lin)
out = lin.result.solution[jls.P1.Isrc.ω]
api    = jls.get_output(sys, lin, lin.result, "P1₊dθ", 1)
manual = out[rm[("P1₊dθ", 1, :Cos)], :] .+ im .* out[rm[("P1₊dθ", 1, :Sin)], :]
@printf("get_output == manual rows: %s\n", isapprox(api, manual))
println("DONE")

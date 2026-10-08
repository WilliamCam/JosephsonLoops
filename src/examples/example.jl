using JosephsonLoops
using ModelingToolkit


loops = [["Lt", "J1", "J2"], ["J2", "Idc"], ["Idc", "R2"]]
squid, _, guesses = build_circuit(process_netlist(loops, ext_flux=[true, false, false]))
squid
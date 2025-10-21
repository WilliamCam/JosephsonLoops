Important: In Julia REPL when the user is asked for an input the initial input is not recognised however the 
           second and all other user input will be recognised.
^--------- Useful to know if using symbolic_assign() or net_edit() in 'sim_model.jl'

******************************************************************************************************************
***	                                           How to use                                                ***
******************************************************************************************************************

1. Run "run_netlist.jl" (Best to run from terminal as Julia REPL is sometimes clunky with readline()) to create 
   a new circuit 
2. If you made any mistakes use the edit function (Accessible thorugh "run_netlist.jl")
3. Run commands from "sim_model.jl" to simulate

******************************************************************************************************************
***	                                      create_netlist.jl                                                ***
******************************************************************************************************************

The file create_netlist.jl is a program used to create and edit circuit 'netlists' which describe circuits in 
terms of loops and components. These 'netlists' are used in sim_model.jl to simulate a circuit.

jdl2 files cannot be edited with a text editor hence to edit a netlist use the edit function in create_netlist.jl

Use run_netlist.jl to use create_netlist.

------------------------------------------------------------------------------------------------------------------
---	                                    Creating a new netlist			       		       ---
------------------------------------------------------------------------------------------------------------------

Total number of loops is limited to 128 loops as numloops is defined as an Int8

..................................................................................................................

CURRENT SOURCES MUST BE ON THE OUTSIDE OF AN EXTERNAL LOOP

When defining loops and components ensure that loops on the right or bottom side of an existing loop have a 
higher loop number than those on the left or top side.

The following loop numbers are good 
  _ _ _   _ _ _   _ _ _ 
|	|	|	|
|   1	|   2	|   3	|
| _ _ _ | _ _ _ | _ _ _ |
|	|	|
|   4	|   5	|
| _ _ _ | _ _ _ |

or
  _ _ _   _ _ _   _ _ _ 
|	|	|	|
|   1	|   3	|   5	|
| _ _ _ | _ _ _ | _ _ _ |
|	|	|
|  2	|  4	|
| _ _ _ | _ _ _ |

These loop numbers are not good and will result in incorrect phase direction and therefore incorrect results
  _ _ _   _ _ _   _ _ _ 
|	|	|	|
|   1	|   4 	|   5	|
| _ _ _ | _ _ _ | _ _ _ |
|	|	|
|   3	|   2	|
| _ _ _ | _ _ _ |

When defining components in a loop separate by line and use the following naming conventions:
Resistor names should start with R				e.g. "R1", "Ra", "Resistor1"
Capacitor names should start with C				e.g. "C1", "Ca", "Capacitor1"
Inductor names should start with L				e.g. "L1", "La", "L_Inductor1"
Josephson Junction names should start with J			e.g. "J1", "Ja", "Josephson1"
Voltage source names should start with V			e.g. "V1", "Va", "V_Source1"
Current source names should start with I			e.g. "I1", "Ia", "I_Source1"

Do not include units or prefixes when defining component parameters, instead enter scientific notation (1e-9)

AC current or voltage source parameters should be entered as 'amplitude, frequency'
	e.g.	What is the amplitude and frequency of V1 (V, Hz)?
		5,60
	^-- This defines an AC Voltage source with an amplitude of 5V and frequency of 60Hz

Entering symbolics for component parameters is supported

..................................................................................................................

Entering symbolics for the external flux vector is not supported

..................................................................................................................
When entering mutually inducted loops enter in the form '1,2,5e-6' where loop 1 and 2 are coupled with mutual 
inductance 5μA/Φ𝜊

------------------------------------------------------------------------------------------------------------------
---                                        Editing an existing netlist				   	       ---
------------------------------------------------------------------------------------------------------------------

It is recommended that a backup file is saved before attempting to edit large netlists as if an error occurs all 
data is corrupted 

..................................................................................................................

Changing component parameters to symbolic values is supported

******************************************************************************************************************
***	                                       model_builder.jl                                                ***
******************************************************************************************************************

The file model_builder.jl is used by sim_model.jl to create a model based on a circuit 'netlist' (which would be 
created using create_netlist.jl)

The user mostly does not directly interact with model_builder.jl as the build function in sim_model.jl automates 
the model building process.

******************************************************************************************************************
***	                                         sim_model.jl                                                  ***
******************************************************************************************************************

The help function in sim_model.jl gives some instructions on how to use the program

Note: using sim_model.jl in terminal is untested and may not function as intended
It is recommended that sim_model.jl is run through a Julia REPL
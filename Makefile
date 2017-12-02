all: julia

julia:
	# run to precompile things that will be needed
	julia genetic.jl 1 < dummy

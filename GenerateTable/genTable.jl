@everywhere push!(LOAD_PATH,"solver")
@everywhere push!(LOAD_PATH,"mdp")
@everywhere using ParallelValueIteration
@everywhere using VerticalCAS

### OPTIONS ###
saveFile = "./VertCAS_qvals_parallel_v7_tauMax40.h5"
###############

# Run in parallel on all available processors
nProcs = nprocs()-1
mdp = VerticalCAS_MDP()
approx_solver = ParallelValueIterationSolver(interp, verbose=true, max_iterations=1,n_procs=nProcs)

# Solve
println("Solving...")
qvals = solve(approx_solver, mdp)

# Write Qvalues to a file
println("Writing Qvalues")
h5open(saveFile, "w") do file
    write(file, "q", qvals) 
end
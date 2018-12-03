module ParallelValueIteration

using Printf
using POMDPs
using POMDPModelTools
using LocalFunctionApproximation
using GridInterpolations
using Distributed
using SharedArrays

import POMDPs: Solver, solve, Policy, action, value 

export
    ValueIterationPolicy,
    ParallelValueIterationSolver,
    solve,
    action,
    value

include("parallelSolver.jl")

end # module
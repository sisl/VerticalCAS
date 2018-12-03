module ParallelValueIteration

using Printf
using POMDPs
using POMDPModelTools
using GridInterpolations
using LocalFunctionApproximation
using Distributed
using SharedArrays
using StaticArrays
using Random

import POMDPs: Solver, solve, Policy, action, value 

export
    ValueIterationPolicy,
    ParallelValueIterationSolver,
    solve,
    action,
    value

include("parallelSolver.jl")

end # module
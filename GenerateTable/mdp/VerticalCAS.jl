module VerticalCAS

using Printf
using POMDPs
using POMDPModelTools
using GridInterpolations
using LocalFunctionApproximation
using StaticArrays

import POMDPs: Solver, solve, Policy, action, value 

include("constants.jl")
include("verticalCAS.jl")
include("transitions.jl")
include("rewards.jl")

end # module
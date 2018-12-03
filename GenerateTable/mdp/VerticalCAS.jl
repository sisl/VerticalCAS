module VerticalCAS

using POMDPs
using POMDPModelTools
using GridInterpolations
using LocalFunctionApproximation
using StaticArrays

include("constants.jl")
include("verticalCAS.jl")
include("transitions.jl")
include("rewards.jl")

end # module
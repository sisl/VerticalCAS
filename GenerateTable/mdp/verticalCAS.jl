export VerticalCAS_MDP

# Define MDP
mutable struct VerticalCAS_MDP <: MDP{stateType, actType}
    hs::Array{Float64,1}
    vowns::Array{Float64,1}
    vints::Array{Float64,1}
    pras::Array{Int64,1}
    resps::Array{Bool,1}
    taus::Array{Float64,1}
    discount_factor::Float64
    accels::Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}
    velRanges::Dict{Int,Tuple{Float64,Float64}}
    allowedTrans::Dict{Int,Vector{Int}}
    currentTau::Float64
end
VerticalCAS_MDP() = VerticalCAS_MDP(hs,vowns,vints,acts,resps,taus,discount_f,accels,velRanges,allowedTrans,0.0)

# Define necessary functions for VerticalCAS MDP
POMDPs.actionindex(::VerticalCAS_MDP, a::actType) = a + 1
POMDPs.actions(mdp::VerticalCAS_MDP)  = acts
POMDPs.discount(mdp::VerticalCAS_MDP) = mdp.discount_factor
#POMDPs.n_actions(::VerticalCAS_MDP)   = length(acts)

function POMDPs.convert_s(::Type{V} where V <: AbstractVector{Float64}, s::stateType, mdp::VerticalCAS_MDP)
    v = [s[1],s[2],s[3],convert(Float64,s[4]),convert(Float64,s[5])] #,s[6]]
    return v
end

function POMDPs.convert_s(::Type{stateType}, v::AbstractVector{Float64}, mdp::VerticalCAS_MDP)
    s = (v[1],v[2],v[3],convert(Int,v[4]), convert(Bool, v[5])) #,v[6])
    return s
end
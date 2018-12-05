export get_belief,get_qval!,Policy,read_policy,evaluate

type Policy
    alpha       :: Matrix{Float64}
    actions     :: Array{Int64,2}
    nactions    :: Int64
    qvals       :: Vector{Float64}

    function Policy(alpha::Matrix{Float64}, actions::Array{Int64,2})
        return new(alpha, actions, size(actions,2), zeros(size(actions,2)))
    end # function Policy
end

function read_policy(actions::Array{Int64,2}, alpha::Matrix{Float64})
    return Policy(alpha, actions)
end # function read_policy

function evaluate(policy::Policy, belief::SparseMatrixCSC{Float64,Int64})
    fill!(policy.qvals, 0.0)
    get_qval!(policy, belief)
    return copy(policy.qvals)
end # function evaluate

function get_qval!(policy::Policy, belief::SparseMatrixCSC{Float64, Int64})
    fill!(policy.qvals, 0.0)
    for iaction in 1:policy.nactions
        for ib in 1:length(belief.rowval)
            policy.qvals[iaction] += belief.nzval[ib] * policy.alpha[belief.rowval[ib], iaction]
        end # for b
    end # for iaction
end # function get_qval!

function get_belief(pstate::Vector{Float64}, grid::RectangleGrid,interp::Bool=false,drl::Bool=false,XandY::Bool=false)
    belief = spzeros(nstates, 1)
    indices, weights = interpolants(grid, pstate)
    if !interp
        largestWeight = 0;
        largestIndex = 0;
        for i = 1:length(weights)
            if weights[i]>largestWeight
                largestWeight = weights[i]
                largestIndex = indices[i]
            end
        end
        indices = largestIndex
        weights = 1.0
    end
    for i = 1:length(indices)
        belief[indices[i]] = weights[i]
    end # for i
    return belief
end # function get_belief
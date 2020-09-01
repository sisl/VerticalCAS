# Defines a solver type
mutable struct ParallelValueIterationSolver{I<:LocalFunctionApproximator} <: Solver
    interp::I
    max_iterations::Int64 # max number of iterations
    belres::Float64 # the Bellman Residual
    verbose::Bool 
    include_Q::Bool
    init_util::Vector{Float64}
    asynchronous::Bool
    n_procs::Int64
end

# Solver with initial value function argument
function ParallelValueIterationSolver(interp::I;
                               max_iterations::Int64 = 100, 
                               belres::Float64 = 1e-3,
                               verbose::Bool = false,
                               include_Q::Bool = true,
                               init_util::Vector{Float64}=Vector{Float64}(undef, 0),
                               asynchronous::Bool = true,
                               n_procs::Int64 = nprocs() - 1)   where
{I<:LocalFunctionApproximator}
    return ParallelValueIterationSolver(interp,max_iterations, belres, verbose, include_Q, init_util, asynchronous, n_procs)
end


# The policy type
mutable struct ValueIterationPolicy <: Policy
    qmat::Matrix{Float64} # Q matrix storing Q(s,a) values
    util::Vector{Float64} # The value function V(s)
    policy::Vector{Int64} # Policy array, maps state index to action index
    action_map::Vector # Maps the action index to the concrete action type
    include_Q::Bool # Flag for including the Q-matrix
    mdp::Union{MDP,POMDP} # uses the model for indexing in the action function
end

# constructor with an optinal initial value function argument
function ValueIterationPolicy(mdp::Union{MDP,POMDP};
                              utility::Vector{Float64}=zeros(n_states(mdp)),
                              policy::Vector{Int64}=zeros(Int64, n_states(mdp)),
                              include_Q::Bool=true)
    ns = n_states(mdp)
    na = length(actions(mdp))
    @assert length(utility) == ns "Input utility dimension mismatch"
    @assert length(policy) == ns "Input policy dimension mismatch"
    action_map = ordered_actions(mdp)
    include_Q ? qmat = zeros(ns,na) : qmat = zeros(0,0)
    return ValueIterationPolicy(qmat, utility, policy, action_map, include_Q, mdp)
end

# constructor for solved q, util and policy
function ValueIterationPolicy(mdp::Union{MDP,POMDP}, q::Matrix{Float64}, util::Vector{Float64}, policy::Vector{Int64})
    action_map = ordered_actions(mdp)
    include_Q = true
    return ValueIterationPolicy(q, util, policy, action_map, include_Q, mdp)
end

# constructor for default Q-matrix
function ValueIterationPolicy(mdp::Union{MDP,POMDP}, q::Matrix{Float64})
    (ns, na) = size(q)
    p = zeros(ns)
    u = zeros(ns)
    for i = 1:ns
        p[i] = argmax(q[i,:])
        u[i] = maximum(q[i,:])
    end
    action_map = ordered_actions(mdp)
    include_Q = true
    return ValueIterationPolicy(q, u, p, action_map, include_Q, mdp)
end


# Function for solving an MDP using parallel value iteration
function solve(solver::ParallelValueIterationSolver,mdp::Union{MDP, POMDP})
    max_iterations = solver.max_iterations
    n_procs = solver.n_procs
    verbose = solver.verbose

    total_time::Float64 = 0.0
    iter_time::Float64 = 0.0

    # Get attributes of interpolator
    # Since the policy object is created by the solver, it directly
    # modifies the value of the interpolator of the created policy
    num_interps::Int = n_interpolating_points(solver.interp)
    interp_points::Vector = get_all_interpolating_points(solver.interp)
    interp_values::Vector = get_all_interpolating_values(solver.interp)

    # Obtain the vector of states by converting the corresponding
    # vector of interpolation points/samples to the state type
    # using the user-provided convert_s function
    println("Computing state space")
    flush(stdout)
    S = statetype(typeof(mdp))
    interp_states = Vector{S}(undef, num_interps)
    for (i,pt) in enumerate(interp_points)
        interp_states[i] = POMDPs.convert_s(S, pt, mdp)
    end

    # init shared utility function and Q-matrix    
    ns = length(interp_states)
    na = length(actions(mdp))

    state_chunks = split_states(ns, n_procs)
    
    println("Sharing with workers")
    flush(stdout)
    
    # intialize the utility and Q-matrix
    if !isempty(solver.init_util)
        @assert length(solver.init_util) == ns "Input utility dimension mismatch"
        init_util = solver.init_util
    else
        init_util = zeros(ns)
    end
    
    include_Q = solver.include_Q
    if include_Q
        init_qmat = zeros(ns, na)
        qmat  = SharedArray{Float64}((ns, na), init = S -> S[localindices(S)] = init_qmat[localindices(S)])
    else
        qmat = nothing 
    end
    
    S = statetype(mdp)
    state_space = interp_states
    shared_states = SharedArray{S}(ns, init = S -> S[localindices(S)] = state_space[localindices(S)])
    util = SharedArray{Float64}(ns, init = S -> S[localindices(S)] = init_util[localindices(S)])    
    maxValue = SharedArray{Float64}(1, init = S -> S[localindices(S)] .= 0.)
    util2 = SharedArray{Float64}(ns, init = S -> S[localindices(S)] = init_util[localindices(S)])
    pool = CachingPool(workers())
    iter_time  = 0.0
    total_time = 0.0
    
    println("Solving chunks")
    flush(stdout)
    
    # Q-matrix to return
    qmats = zeros(ns*length(mdp.taus),na)
    
    for i =1:length(mdp.taus)
        iter_time = @elapsed begin
            maxValue[1] = 0.
            copyto!(util2,util)
            results = pmap(x -> solve_chunk(mdp, mdp.taus[i], solver.interp,shared_states, util, util2, qmat, include_Q, maxValue, x), pool, state_chunks)
        end
        total_time += iter_time
        verbose ? @printf("[Tau: %-4d] Max Absolute Q-Value: %10.5f | iteration runtime: %10.3f s, (%10.3G s total)\n", mdp.taus[i], maxValue[1], iter_time, total_time) : nothing
        verbose ? flush(stdout) : nothing
        
        qmats[(1+(i-1)*ns):(i*ns),:] = deepcopy(convert(Array{Float64,2},qmat))
    end # main iteration loop
        
    return qmats
end

# Updates shared utility and Q-Matrix using gauss-seidel value iteration (asynchronous)
function solve_chunk(mdp::M, 
                    tau::Float64,
                    interp::LocalFunctionApproximator,
                    states::SharedArray{S, 1}, 
                    util::SharedArray{Float64, 1}, 
                    util2::SharedArray{Float64,1},
                    qmat::SharedArray{Float64, 2}, 
                    include_Q::Bool,
                    maxValue::SharedArray{Float64, 1},
                    state_indices::UnitRange{Int64}
                    ) where {M <: Union{MDP, POMDP}, S}

    discount_factor = discount(mdp)
    mdp.currentTau = tau
    for istate=state_indices
        s = states[istate]
        sub_aspace = actions(mdp, s)
        if isterminal(mdp, s)
            util[istate] = 0.0
        else
            old_util = util[istate] # for maxValue
            max_util = -Inf
            for a in sub_aspace
                iaction = actionindex(mdp, a)
                dist = transition(mdp, s, a) # creates distribution over neighbors
                u = 0.0
                r = reward(mdp, s, a)
                for (sp, p) in weighted_iterator(dist)
                    p == 0.0 ? continue : nothing # skip if zero prob                    
                    u += p*r
                    # Only interpolate sp if it is non-terminal
                    if !isterminal(mdp,sp)
                        sp_point = POMDPs.convert_s(Vector{Float64}, sp, mdp)
                        u += p*discount_factor*interpolate(interp.grid,util2,sp_point)
                    end
                end
                new_util = u
                if new_util > max_util
                    max_util = new_util
                end
                include_Q ? (qmat[istate, iaction] = new_util) : nothing
            end # action
            
            
            util[istate] = max_util
            diff = abs(max_util)
            
            diff > maxValue[1] ? (maxValue[1] = diff) : nothing
        end
    end
    return 
end

# Split the states into separate chunks based on number of processors
function split_states(ns::Int64, n_procs::Int64)
    state_chunks = Vector{UnitRange{Int64}}(undef, n_procs)
    stride = div(ns, n_procs)
    for j=0:n_procs-1
        sj = j*stride + 1
        ej = (j + 1)*stride
        if j == n_procs-1
            ej = ns
        end
        state_chunks[j+1] = sj:ej
    end 
    return state_chunks
end
    
function action(policy::ValueIterationPolicy, s::S) where S
    sidx = stateindex(policy.mdp, s)
    aidx = policy.policy[sidx]
    return policy.action_map[aidx]
end

function value(policy::ValueIterationPolicy, s::S) where S
    sidx = stateindex(policy.mdp, s)
    policy.util[sidx]
end

@everywhere using POMDPs
@everywhere using POMDPModelTools
@everywhere using GridInterpolations
@everywhere using LocalFunctionApproximation

# Wrote parallel version of LocalApproximationValueIteration, no need to include official version
#using LocalApproximationValueIteration

@everywhere using Printf
@everywhere using Distributed
@everywhere using SharedArrays
@everywhere using Random
@everywhere using StaticArrays
@everywhere using HDF5
@everywhere import POMDPs: Solver, solve, Policy, action, value

### OPTIONS ###
saveFile = "VertCAS_qvals_parallel_v7_tauMax40.h5"
nProcs = 7
###############

############ SOLVER RELATED FUNCTIONS ####################################
@everywhere mutable struct ParallelLocalApproximationValueIterationSolver{I<:LocalFunctionApproximator} <: Solver
    interp::I
    max_iterations::Int64 # max number of iterations
    belres::Float64 # the Bellman Residual
    verbose::Bool 
    include_Q::Bool
    init_util::Vector{Float64}
    asynchronous::Bool
    n_procs::Int64
end

@everywhere function ParallelLocalApproximationValueIterationSolver(interp::I;
                               max_iterations::Int64 = 100, 
                               belres::Float64 = 1e-3,
                               verbose::Bool = false,
                               include_Q::Bool = true,
                               init_util::Vector{Float64}=Vector{Float64}(undef, 0),
                               asynchronous::Bool = true,
                               n_procs::Int64 = nprocs() - 1)   where
{I<:LocalFunctionApproximator}
    return ParallelLocalApproximationValueIterationSolver(interp,max_iterations, belres, verbose, include_Q, init_util, asynchronous, n_procs)
end


# The policy type
@everywhere mutable struct ValueIterationPolicy <: Policy
    qmat::Matrix{Float64} # Q matrix storing Q(s,a) values
    util::Vector{Float64} # The value function V(s)
    policy::Vector{Int64} # Policy array, maps state index to action index
    action_map::Vector # Maps the action index to the concrete action type
    include_Q::Bool # Flag for including the Q-matrix
    mdp::Union{MDP,POMDP} # uses the model for indexing in the action function
end

# constructor with an optinal initial value function argument
@everywhere function ValueIterationPolicy(mdp::Union{MDP,POMDP};
                              utility::Vector{Float64}=zeros(n_states(mdp)),
                              policy::Vector{Int64}=zeros(Int64, n_states(mdp)),
                              include_Q::Bool=true)
    ns = n_states(mdp)
    na = n_actions(mdp)
    @assert length(utility) == ns "Input utility dimension mismatch"
    @assert length(policy) == ns "Input policy dimension mismatch"
    action_map = ordered_actions(mdp)
    include_Q ? qmat = zeros(ns,na) : qmat = zeros(0,0)
    return ValueIterationPolicy(qmat, utility, policy, action_map, include_Q, mdp)
end

# constructor for solved q, util and policy
@everywhere function ValueIterationPolicy(mdp::Union{MDP,POMDP}, q::Matrix{Float64}, util::Vector{Float64}, policy::Vector{Int64})
    action_map = ordered_actions(mdp)
    include_Q = true
    return ValueIterationPolicy(q, util, policy, action_map, include_Q, mdp)
end

# constructor for default Q-matrix
@everywhere function ValueIterationPolicy(mdp::Union{MDP,POMDP}, q::Matrix{Float64})
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


@everywhere function solve(solver::ParallelLocalApproximationValueIterationSolver,mdp::Union{MDP, POMDP})
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
    ns = length(interp_states)#n_states(mdp)
    na = n_actions(mdp)

    #chunk_batches = split_states(ns, n_procs,n_batches)
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

# updates shared utility and Q-Matrix using gauss-seidel value iteration (asynchronous)
@everywhere function solve_chunk(mdp::M, 
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
            for a in [0,1,2,3,4,5,6,7,8]#sub_aspace
                iaction = actionindex(mdp, a)
                dist = transition(mdp, s, a) # creates distribution over neighbors
                u = 0.0
                for (sp, p) in weighted_iterator(dist)
                    p == 0.0 ? continue : nothing # skip if zero prob
                    r = reward(mdp, s, a, sp)
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
            diff = abs(max_util) # - old_util)
            
            diff > maxValue[1] ? (maxValue[1] = diff) : nothing
        end
    end # state loop
    return 
end

@everywhere function split_states(ns::Int64, n_procs::Int64)
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
    
@everywhere function action(policy::ValueIterationPolicy, s::S) where S
    sidx = stateindex(policy.mdp, s)
    aidx = policy.policy[sidx]
    return policy.action_map[aidx]
end

@everywhere function value(policy::ValueIterationPolicy, s::S) where S
    sidx = stateindex(policy.mdp, s)
    policy.util[sidx]
end
##################################################################################


println("Defining Problem...")
flush(stdout)
# Advisory indices
@everywhere begin
    
    # Define constants
    COC=0
    DNC=1
    DND=2
    DES1500 = 3
    CL1500 = 4
    SDES1500=5
    SCL1500=6
    SDES2500=7
    SCL2500=8

    # State Type:
    stateType = Tuple{Float64,Float64,Float64,Int,Bool}
    actType = Int

    # Default parameters
    discount_f = 1.0
    hMin = -8000.0
    hMax = 8000.0
    numH = 641

    vMin = -100.
    vMax = 100.
    numV = 51

    tauMin = 0.0
    tauMax = 40.0
    numTau = 41
    acts = [0,1,2,3,4,5,6,7,8]
    resps = [false,true]

    vels = vcat(LinRange(-100,-60,5),LinRange(-50,-35,4),LinRange(-30,30,21),LinRange(35,50,4),LinRange(60,100,5))
    hs   = vcat(LinRange(-8000,-4000,5),LinRange(-3000,-1250,8),LinRange(-1000,-800,3),LinRange(-700,-150,12),LinRange(-100,100,9),LinRange(150,700,12),LinRange(800,1000,3),LinRange(1250,3000,8),LinRange(4000,8000,5))
    vowns = vels #LinRange(vMin,vMax,6)
    vints = vels #LinRange(vMin,vMax,6)
    taus  = LinRange(tauMin,tauMax,numTau)

    grid = RectangleGrid(hs,vowns,vints,acts,resps) # Create the interpolating grid
    interp = LocalGIFunctionApproximator(grid)  # Create the local function approximator using the grid

    accels = Dict(COC=>([0.5,0.25,0.25],[0.0,3.0,-3.0]),
                      DNC=>([0.5,0.25,0.25],[-8.33,-9.33,-7.33]),
                      DND=>([0.5,0.25,0.25],[8.33,9.33,7.33]),
                      DES1500=>([0.5,0.25,0.25],[-8.33,-9.33,-7.33]),
                      CL1500=>([0.5,0.25,0.25],[8.33,9.33,7.33]),
                      SDES1500=>([0.5,0.25,0.25],[-10.7,-11.7,-9.7]),
                      SCL1500=>([0.5,0.25,0.25],[10.7,11.7,9.7]),
                      SDES2500=>([0.5,0.25,0.25],[-10.7,-11.7,-9.7]),
                      SCL2500=>([0.5,0.25,0.25],[10.7,11.7,9.7]))

    velRanges = Dict(COC=>(-100.0,100.0),
                    DNC=>(0.0,100.0),
                    DND=>(-100.0,0.0),
                    DES1500=>(-25.0,100.0),
                    CL1500=>(-100.0,25.0),
                    SDES1500=>(-25.0,100.0),
                    SCL1500=>(-100.0,25.0),
                    SDES2500=>(-41.67,100.0),
                    SCL2500=>(-100.0,41.67))
    
    allowedTrans = Dict(COC=>[1,1,1,1,1,0,0,0,0],
                           DNC=>[1,1,1,1,1,0,0,0,0],
                           DND=>[1,1,1,1,1,0,0,0,0],
                           DES1500=>[1,1,1,1,1,1,1,0,0],
                           CL1500=>[1,1,1,1,1,1,1,0,0],
                           SDES1500=>[1,1,1,1,1,1,1,1,1],
                           SCL1500=>[1,1,1,1,1,1,1,1,1],
                           SDES2500=>[1,1,1,1,1,1,1,1,1],
                           SCL2500=>[1,1,1,1,1,1,1,1,1])
    
    # Define MDP
    mutable struct VerticalCAS <: MDP{stateType, actType}
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
    VerticalCAS() = VerticalCAS(hs,vowns,vints,acts,resps,taus,discount_f,accels,velRanges,allowedTrans,0.0)

    # Define necessary functions for VerticalCAS MDP
    POMDPs.actionindex(::VerticalCAS, a::actType) = a + 1
    POMDPs.actions(mdp::VerticalCAS)  = acts
    POMDPs.discount(mdp::VerticalCAS) = mdp.discount_factor
    POMDPs.n_actions(::VerticalCAS)   = length(acts)

    function POMDPs.convert_s(::Type{V} where V <: AbstractVector{Float64}, s::stateType, mdp::VerticalCAS)
        v = [s[1],s[2],s[3],convert(Float64,s[4]),convert(Float64,s[5])] #,s[6]]
        return v
    end

    function POMDPs.convert_s(::Type{stateType}, v::AbstractVector{Float64}, mdp::VerticalCAS)
        s = (v[1],v[2],v[3],convert(Int,v[4]), convert(Bool, v[5])) #,v[6])
        return s
    end

    function POMDPs.transition(mdp::VerticalCAS, s::stateType, ra::actType)
        h = s[1]; vown = s[2]; vint = s[3]; pra = s[4]; prev_resp = s[5]#; tau = s[6]
        nextStates = MVector{18, stateType}(undef)
        nextProbs = @MVector(zeros(18))
        ind=1

        prob_resp = 0.0
        if (ra==COC) .| (prev_resp .& (pra==ra))
            prob_resp=1.0
        elseif (pra==COC) #.& (ra!=COC)
            prob_resp = 1.0/(1+5.)
        elseif sameSense(pra,ra)
            prob_resp = 1.0/(1+3.)
        else
            prob_resp = 1.0/(1+5.)
        end

        next_resp = true
        next_pra  = ra
        ownProbs, ownAccels = mdp.accels[ra]
        intProbs, intAccels = mdp.accels[COC]
        for i = 1:3
            for j = 1:3
                next_h,next_vown,next_vint = dynamics(h,vown,vint,ownAccels[i],intAccels[j],ra,mdp)
                nextStates[ind] = (next_h,next_vown,next_vint,next_pra,next_resp)#,next_tau)
                nextProbs[ind]  = prob_resp*ownProbs[i]*intProbs[j]
                ind+=1
            end
        end

        if prob_resp<1.0
            prob_notResp = 1.0-prob_resp
            next_resp=false
            ownProbs, ownAccels = mdp.accels[COC]
            for i = 1:3
                for j = 1:3
                    next_h,next_vown,next_vint = dynamics(h,vown,vint,ownAccels[i],intAccels[j],COC,mdp)
                    nextStates[ind] = (next_h,next_vown,next_vint,next_pra,next_resp)#,next_tau)
                    nextProbs[ind]  = prob_notResp*ownProbs[i]*intProbs[j]
                    ind+=1
                end
            end 
        end

        return SparseCat(nextStates,nextProbs)
    end

    function dynamics(h::Float64,vown::Float64,vint::Float64,ownAccel::Float64, intAccel::Float64, ra::Int, mdp::VerticalCAS)
        vLow, vHigh = mdp.velRanges[ra]
        if (vLow >= vown) .| (vHigh <= vown)
            ownAccel = 0
        elseif vLow > vown + ownAccel
            ownAccel = vLow-vown
        elseif vHigh < vown + ownAccel
            ownAccel = vHigh-vown
        end

        vLow, vHigh = mdp.velRanges[COC]
        if (vLow >= vint) .| (vHigh <= vint)
            intAccel = 0
        elseif vLow > vint + intAccel
            intAccel = vLow-vint
        elseif vHigh < vint + intAccel
            intAccel = vHigh-vint
        end

        next_h = h-vown-0.5*ownAccel+vint+0.5*intAccel
        next_vown = vown+ownAccel
        next_vint = vint+intAccel
        return next_h,next_vown,next_vint
    end

    sameSense(pra::Int, ra::Int) = mod(pra,2)==mod(ra,2)
    downSense(ra::Int) = (ra>0) .& (mod(ra,2)==1)
    upSense(ra::Int) = (ra>0) .& (mod(ra,2)==0)

    function POMDPs.reward(mdp::VerticalCAS, s::stateType, ra::actType)
        h = s[1]; vown = s[2]; vint = s[3]; pra = s[4]; resp = s[5]; #tau = s[6]
        tau = mdp.currentTau
        r = 0.0
        sep = abs(h)
        closure = abs(vint-vown)
        crossing = ((h<0) .& downSense(ra)) .| ((h>0) .& upSense(ra))
        deltaVown = 0.0
        corrective = false
        preventitive = false
        weakening = false
        strengthening = false
        reversal = false
        
        if ra>COC
            if pra>COC
                reversal = mod(pra,2)!=mod(ra,2)
            end
            if !reversal
                weakening = pra>ra
                strengthening = pra<ra
            end
            vLow, vHigh = mdp.velRanges[ra]
            corrective = (vown>vLow) .& (vown < vHigh)
            preventitive = !corrective
            if corrective
                if downSense(ra)
                    deltaVown = abs(vLow-vown)
                else
                    deltaVown = abs(vHigh-vown)
                end
            end
        end
        lolo = (ra==DNC) .| (ra==DND)

        if (sep<=175) .& (tau==0)
            r-=1.0
        end
    
    
        if mdp.allowedTrans[pra][ra+1]==0
            r-=1.0
        end
    
        if crossing
            if preventitive
                r-=1.0
            end
            if sep>500
                r-=0.01
            end
        end
        if corrective
            r-=1e-5
            if (sep>650) .& (closure<2000.0/60.0)
                r-= 0.1
            end
            if (sep>1000) .& (closure<4000.0/60.0)
                r-=0.03
            end
        elseif preventitive
            if (sep>650) .& (closure<2000.0/60.0)
                r-=0.01
            end
        end
        if reversal
            r-= 8e-3
        end
        if strengthening
            r-=5e-3
        end
        if weakening
            r-=1e-3
        end
        if lolo
            r-=1e-4
            if closure > 3000.0/60.0
                r-=5e-4
            end
        elseif (ra!=COC) .& (closure > 3000.0/60.0)
            r-=1.5e-3
        end
        if closure < 3000.0/60.0
            r-=2.3e-3
        end
        if ra==COC
            r+=1e-9
        else
            r-=3e-5*deltaVown
            if closure > 3000.0/60.0
                r-=1.5e-3
            end
        end
        return r
    end
end

println("Making MDP")
flush(stdout)

mdp = VerticalCAS()
approx_solver = ParallelLocalApproximationValueIterationSolver(interp, verbose=true, max_iterations=1,n_procs=nProcs)

println("Solving...")
flush(stdout)

qvals = solve(approx_solver, mdp)

println("Writing Qvalues")
flush(stdout)
h5open(saveFile, "w") do file
    write(file, "q", qvals) 
end
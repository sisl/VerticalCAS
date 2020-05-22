@everywhere push!(LOAD_PATH,"mdp")
@everywhere using LocalFunctionApproximation
@everywhere using GridInterpolations
@everywhere using POMDPs
@everywhere using SparseArrays
@everywhere using VerticalCAS
@everywhere using HDF5
@everywhere using POMDPModelTools
@everywhere using Printf


### OPTIONS ###
saveFile = "/raid/kjulian3/VertCAS/Qtables/VertCAS_noResp_newTrans_v5.h5"
nTau0=10   # Number of seconds at tau=0
maxTau=40  # Max tau value
###############

@everywhere function compute_helper(states,n_states,mdps,grid,a)
    # Initialize the row, column and z vectors to be large vectors.
    # We truncate these vectors at the end when we return the sparse array
    rval = zeros(Int32,n_states*100)
    cval = zeros(Int32,n_states*100)
    zval = zeros(Float32,n_states*100)
    
    # Initialize reward vector, 1 for each MDP, since the reward is different for 
    # different tau values
    rews = []
    for t = 1:length(mdps)
        push!(rews,zeros(n_states))
    end
    
    # Printing
    index=1
    fracBase=0.2
    frac = fracBase
    
    ## Iterate through all states
    for (i,s) in enumerate(states)
        if i/n_states>=frac
            print(round(frac*100))
            println("% Complete")
            frac+=fracBase
        end

        ## Compute rewards
        for j=1:length(mdps)
            rews[j][i] = reward(mdps[j],s,a)
        end

        ## Compute transitions, just use the first MDP (dynamics are same for all MDPs!)
        dist = transition(mdps[1], s, a)
        for (sp, p) in weighted_iterator(dist)
            if p >0.0
                sp_point = POMDPs.convert_s(Vector{Float64}, sp, mdps[1])
                sps, probs = interpolants(grid, sp_point)
                for (spi, probi) in zip(sps,probs)
                    rval[index] = i
                    cval[index] = spi
                    zval[index] = probi*p
                    index += 1
                end
            end
        end
    end
    
    # Create sparse transition matrix, removing unused space in row, column and z vectors
    trans = sparse(rval[1:index-1],cval[1:index-1],zval[1:index-1],n_states,n_states)
    return (trans, rews)
end

# Bellman Update
@everywhere function compute_Qa(r,gam,trans,U)
    return r + gam*trans*U
end

# Calls compute_helper to compute the transition sparse array and reward vectors
# Compute the arrays for each action in parallel using processors as available
function compute_trans_reward(mdps::Array{Union{MDP, POMDP},1},interp::LocalFunctionApproximator)
    ## Dictionaries to populate
    t = Dict()
    rews = Dict()
    rc = Dict()
    
    ## Compute states
    n_states = n_interpolating_points(interp)
    interp_points = get_all_interpolating_points(interp)
    S = statetype(typeof(mdps[1]))
    interp_states = Vector{S}(undef, n_states)
    for (i,pt) in enumerate(interp_points)
        interp_states[i] = POMDPs.convert_s(S, pt, mdps[1])
    end
    
    ## Loop through all actions in parallel
    # Split the tasks between workers
    for (ai,a) in enumerate(actions(mdps[1]))
        rc[ai] = remotecall(compute_helper, mod(ai,nprocs()-1)+2,interp_states,n_states,mdps,interp.grid,a)
    end
    # Call the workz
    for (ai,a) in enumerate(actions(mdps[1]))
        (t[ai], rews[ai]) = fetch(rc[ai])
    end
    return (t, rews)
end

function computeQ(mdps::Array{Union{MDP, POMDP},1},interp,nTau0)
    (trans, rews) = compute_trans_reward(mdps,interp);

    # Initialize Q and U vectors
    ns = length(rews[1][1])
    na = length(actions(mdps[1]))
    nt = length(mdps)
    taus = 0:(nt-1)
    gam = discount(mdps[1])
    U = zeros(ns)
    Q = zeros(ns,na)
    Q_out = zeros(ns*nt,na)
    Q_rc = Dict()

    ## Warm start for U at tau=0
    for i=1:nTau0
        @printf("Warm up: %d/%d\n",i,nTau0)
        for ai = 1:na
            Q_rc[ai] = remotecall(compute_Qa,mod(ai,nprocs()-1)+2,rews[ai][1],gam,trans[ai],U)
        end
        for ai = 1:na
            Q[:,ai] = fetch(Q_rc[ai])
        end
        U = maximum(Q,dims=2)
    end

    for i=1:nt
        @printf("Tau %d\n",taus[i])
        for ai = 1:na
            Q_rc[ai] = remotecall(compute_Qa,mod(ai,nprocs()-1)+2,rews[ai][i],gam,trans[ai],U)
        end
        for ai = 1:na
            Q[:,ai] = fetch(Q_rc[ai])
        end
        U = maximum(Q,dims=2)
        Q_out[1+(i-1)*ns:(i*ns),:] = deepcopy(Q)
    end
    return Q_out
end

mdps = Array{Union{MDP, POMDP},1}()
for tau = 0:maxTau
    push!(mdps, VerticalCAS_MDP())
    mdps[end].currentTau = tau
end
@time Q_out = computeQ(mdps,interp,nTau0)


println("Writing Qvalues")
h5open(saveFile, "w") do file
    write(file, "q", Q_out) 
end

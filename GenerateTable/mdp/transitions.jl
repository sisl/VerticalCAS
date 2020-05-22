# State transition function
function POMDPs.transition(mdp::VerticalCAS_MDP, s::stateType, ra::actType)
    h = s[1]; vown = s[2]; vint = s[3]; pra = s[4];
    
    # Computation is faster when using vector of static size
    nextStates = MVector{9, stateType}(undef)
    nextProbs = @MVector(zeros(18))
    next_pra = ra
    ind=1

    # Compute probabilities of next states using sigma point sampling
    ownProbs, ownAccels = mdp.accels[pra]
    intProbs, intAccels = mdp.accels[-1]
    for i = 1:3
        for j = 1:3
            next_h,next_vown,next_vint = dynamics(h,vown,vint,ownAccels[i],intAccels[j],pra,mdp)
            nextStates[ind] = (next_h,next_vown,next_vint,next_pra)
            nextProbs[ind]  = ownProbs[i]*intProbs[j]
            ind+=1
        end
    end

    return SparseCat(nextStates,nextProbs)
end

# Dynamic equations
function dynamics(h::Float64,vown::Float64,vint::Float64,ownAccel::Float64, intAccel::Float64, ra::Int, mdp::VerticalCAS_MDP)
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
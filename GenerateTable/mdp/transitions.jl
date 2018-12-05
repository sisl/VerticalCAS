# State transition function
function POMDPs.transition(mdp::VerticalCAS_MDP, s::stateType, ra::actType)
    h = s[1]; vown = s[2]; vint = s[3]; pra = s[4]; prev_resp = s[5]
    
    # Computation is faster when using vector of static size
    nextStates = MVector{18, stateType}(undef)
    nextProbs = @MVector(zeros(18))
    ind=1

    # Compute probability that new advisory will be followed
    prob_resp = 0.0
    if (ra==COC) .| (prev_resp .& (pra==ra))
        prob_resp=1.0
    elseif (pra==COC)
        prob_resp = 1.0/(1+5.)
    elseif sameSense(pra,ra)
        prob_resp = 1.0/(1+3.)
    else
        prob_resp = 1.0/(1+5.)
    end

    # Compute probabilities of next states using sigma point sampling
    next_resp = true
    next_pra  = ra
    ownProbs, ownAccels = mdp.accels[ra]
    intProbs, intAccels = mdp.accels[COC]
    for i = 1:3
        for j = 1:3
            next_h,next_vown,next_vint = dynamics(h,vown,vint,ownAccels[i],intAccels[j],ra,mdp)
            nextStates[ind] = (next_h,next_vown,next_vint,next_pra,next_resp)
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
                nextStates[ind] = (next_h,next_vown,next_vint,next_pra,next_resp)
                nextProbs[ind]  = prob_notResp*ownProbs[i]*intProbs[j]
                ind+=1
            end
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
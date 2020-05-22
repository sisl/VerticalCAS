# Helper functions
sameSense(pra::Int, ra::Int) = mod(pra,2)==mod(ra,2)
downSense(ra::Int) = (ra>0) .& (mod(ra,2)==1)
upSense(ra::Int) = (ra>0) .& (mod(ra,2)==0)

# Reward function for VerticalCAS MDP
function POMDPs.reward(mdp::VerticalCAS_MDP, s::stateType, ra::actType)
    h = s[1]; vown = s[2]; vint = s[3]; pra = s[4]
    tau = mdp.currentTau
    r = 0.0
    sep = abs(h)
    sepAtTau0 = abs(h+vint*tau - vown*tau)
    
    # Set closure to 0 if aircraft are diverging!
    closure = vint - vown
    if (closure>0) .& (h>0)
        closure = 0
    elseif (closure<0) .& (h<0)
        closure = 0
    else
        closure = abs(closure)
    end
    
    crossing = ((h<-100) .& downSense(ra)) .| ((h>100) .& upSense(ra))
    deltaVown = 0.0
    corrective = false
    preventive = false
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
        corrective = (vown>vLow) .& (vown < vHigh) .& (ra != COC)
        preventive = !corrective .& (ra != COC)
        if corrective
            if downSense(ra)
                deltaVown = abs(vLow-vown)
            else
                deltaVown = abs(vHigh-vown)
            end
        end
    end
    # Level-off, Level-off
    lolo = ((ra==DNC) .| (ra==DND)) .& corrective
    # Maintain Vertical Speed
    mvs = (ra!=COC) .& preventive

    
    # Collision penalty
    r -= 1.5*exp(-sep/200.0)*exp(-tau/3.0)
    
    # Penalty for ending advisory too early
    if (ra==COC) .& (pra!=COC)
        r -= 1.0*exp(-sepAtTau0/100.0)*exp(-tau/10.0)
    elseif (pra==COC) .& ((ra==COC) .| (ra==DNC) .| (ra==DND))
        # Penalty for COC when a collision is close approaching
        r -= 0.5*exp(-sepAtTau0/200.0)*exp(-tau/5.0)
    end
    
    #if (sep<=175) .& (tau==0)
    #    r-=1.0*2.0 # collision penalty
    #end
    
    #if (sep<=100) .& (tau<10) #.& (ra==COC)
    #    r-=0.1
    #end
        
    
    if mdp.allowedTrans[pra][ra+1]==0
        r-=1.0*10.0 # was 2
    end
    
    
    if crossing
        if preventive .& (sep>650)
            r-=1 #Remove this penalty?
        end
        if sep>500
            r-=0.01
        end
        if abs(vown)>500.0/60.0
            if ((sign(vown)==1) .& (downSense(ra))) .| ((sign(vown)==-1) .& (upSense(ra)))
                r -= 4e-4*deltaVown
            end
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
    elseif preventive
        if (sep>650) .& (closure<2000.0/60.0)
            r-=0.01
        end
    end
    
    if reversal
        r-= 8e-3 # Reversal penalty
        
        #### NEW! Penalize reversals in NMAC region - with pilot delay, better to be consistent
        #if sep<100
        #    r-=0.25
        #end
    end
    if strengthening
        r-=5e-3
    end
    if weakening
        r-=1e-3
    end
    
    if lolo .| mvs
        r-=1.0e-4
        if (closure > 3000.0/60.0)
            r-=5e-4
        end
    elseif (ra!=COC) .& (closure > 3000.0/60.0)
        r-=1.5e-3
    end
    if (closure < 3000.0/60.0) .& (ra!=COC)
        r-=2.3e-5
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
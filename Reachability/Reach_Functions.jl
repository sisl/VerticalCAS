function getNextSetBounds!(cellBounds,accelBounds)
    #= Get an over-approximated region of space reachable from each cell
    Inputs: 
        hBounds (array of floats): lower/upper bounds on h for each cell
        vownBounds (array of floats): lower/upper bounds on vown for each cell
        vintBounds (array of floats): lower/upper bounds on vint for each cell
        a0min (array of floats): lower bounds on a0 (ownship acceleration) for each cell
        a0max (array of floats): Upper bounds on a0 (ownship acceleration) for each cell
        a1min (array of floats): lower bounds on a1 (intruder acceleration) for each cell
        a1max (array of floats): Upper bounds on a1 (intruder acceleration) for each cell
    Outputs:
        None (hBounds, vintBounds, and vownBounds are changed to be the boundaries of reachable space for each cell
    =#
    
    # Unpack values
    hBounds,vownBounds,vintBounds = cellBounds
    a0min,a0max,a1min,a1max,v0min,v0max = accelBounds
    
    indsMinLow = vownBounds[:,1] .< v0min
    indsMinHi = vownBounds[:,1] .> v0max
    indsMaxLow = vownBounds[:,2] .< v0min
    indsMaxHi = vownBounds[:,2] .> v0max
    
    
    tAccelMin = zeros(Float32,size(hBounds)[1])
    tAccelMax = zeros(Float32,size(hBounds)[1])
    if v0min < v0max
        tAccelMin[indsMinLow] = (v0min .- vownBounds[indsMinLow,1])./a0min
        tAccelMin[indsMinHi] = (v0max .- vownBounds[indsMinHi,1])./a0min
        tAccelMin[tAccelMin.>1] .= 1
        
        tAccelMax[indsMaxLow] = (v0min .- vownBounds[indsMaxLow,2])./a0max
        tAccelMax[indsMaxHi] = (v0max .- vownBounds[indsMaxHi,2])./a0max
        tAccelMax[tAccelMax.>1] .= 1
    else
        tAccelMax.=1
        tAccelMin.=1
    end
    
    hBounds[:,1] += vintBounds[:,1] .+ 0.5*a1min - vownBounds[:,2] - (0.5*a0max.*tAccelMin)
    hBounds[:,2] += vintBounds[:,2] .+ 0.5*a1max - vownBounds[:,1] - (0.5*a0min.*tAccelMax)
    vownBounds[:,1] += a0min.*tAccelMin
    vownBounds[:,2] += a0max.*tAccelMax
    vintBounds[:,1] .+= a1min
    vintBounds[:,2] .+= a1max
end
    
function getCutpoints(XS)
    xStepsAll = XS[2:end]-XS[1:end-1]
    xSteps = Array{Float32,1}()
    xNums = Array{Int32,1}()
    xCuts = [XS[1]]
    currentStep = -1
    for (i,step) in enumerate(xStepsAll)
        if step != currentStep
            currentStep = step
            xSteps = vcat(xSteps,[step])
            xNums = vcat(xNums,[i])
            if length(xNums)>1
                xCuts = vcat(xCuts,[xCuts[end]+(xNums[end]-xNums[end-1])*xSteps[end-1]])
            end
        end
    end
    return (xSteps, xCuts, xNums)
end


function getIndexBounds(cellBounds)
    #= Compute which cells overlap with region of state space for each region
    Inputs:
        hBoundsNext (array of floats): Lower/upper bounds on h for each region
        vownBoundsNext (array of floats): Lower/upper bounds on vown for each region
        vintBoundsNext (array of floats): Lower/upper bounds on vint for each region
    Outputs:
        hInds (array of ints): First/last index of h-array that overlap with the region
        vownInds (array of ints): First/last index of vown-array that overlap with the region
        vintInds (array of ints): First/last index of vint-array that overlap with the region
    =#
    
    # Only use cell if it is inside, not just touching the boundary. This small adjustment 
    # prevents cells from being added that just touch the boundary of the region
        
    hBoundsNext, vownBoundsNext, vintBoundsNext = cellBounds
    vintBoundsNext[:,1].+=1e-6
    vintBoundsNext[:,2].-=1e-6
    hBoundsNext[:,1].+=1e-6
    hBoundsNext[:,2].-=1e-6
    vownBoundsNext[:,1].+=1e-6
    vownBoundsNext[:,2].-=1e-6
        
    # The X/Y arrays of the grid are not uniform, but because it is piecewise linearly spaced, we
    # can compute the indices more efficiently than a brute force search
    hInds = ones(Int32,size(hBoundsNext)[1],2)
    hMap = hInds.>0
    hSteps, hCuts, hNums = getCutpoints(HS)
    for i=1:length(hSteps)
        hMap = hMap .& (hBoundsNext.>hCuts[i])
        hInds[hMap] = Int32.(div.(hBoundsNext[hMap].+(hSteps[i]*hNums[i]-hCuts[i]),hSteps[i]))
    end
    
    vownInds = ones(Int32,size(vownBoundsNext)[1],2)
    vownMap = vownInds.>0
    vownSteps, vownCuts, vownNums = getCutpoints(VOWNS)
    for i=1:length(vownSteps)
        vownMap = vownMap .& (vownBoundsNext.>vownCuts[i])
        vownInds[vownMap] = Int32.(div.(vownBoundsNext[vownMap].+(vownSteps[i]*vownNums[i]-vownCuts[i]),vownSteps[i]))
    end
    
    vintInds = ones(Int32,size(vintBoundsNext)[1],2)
    vintMap = vintInds.>0
    vintSteps, vintCuts, vintNums = getCutpoints(VINTS)
    for i=1:length(vintSteps)
        vintMap = vintMap .& (vintBoundsNext.>vintCuts[i])
        vintInds[vintMap] = Int32.(div.(vintBoundsNext[vintMap].+(vintSteps[i]*vintNums[i]-vintCuts[i]),vintSteps[i]))
    end
    
    hInds[hInds.>NUMH] .= NUMH
    vownInds[vownInds.>NUMVOWN] .= NUMVOWN
    vintInds[vintInds.>NUMVINT] .= NUMVINT

    return (hInds,vownInds,vintInds)
end


function getReachDynamics(cellBoundInds)
    #= Convert the first/last indices for each dimension to a sparse
       representation of the reachable dynamics. Output r has length the same as the number of cells.
       The cells reachable from cell i are given by c[r[i]:r[i+1]-1]
    
    This function is slow due to for-loops
    
    Inputs:
        xBoundInds (array of ints): First/last index of x-array that overlap with the region
        yBoundInds (array of ints): First/last index of y-array that overlap with the region
        psiBoundInds (array of ints): First/last index of psi-array that overlap with the region
    Outputs:
        r (array of ints): Pointers to the c array to determine where to begin reading
        c (array of ints): Array of cell indices
    =#
    
    hBoundInds,vownBoundInds,vintBoundInds = cellBoundInds
    
    dH = Int32(maximum(hBoundInds[:,2]-hBoundInds[:,1]))
    dVown = Int32(maximum(vownBoundInds[:,2]-vownBoundInds[:,1]))
    dVint = Int32(maximum(vintBoundInds[:,2]-vintBoundInds[:,1]))
    
    r = Vector{Int32}()
    c = Vector{Int32}()
    nr = size(hBoundInds)[1]
    inds = Int32(1):nr
    
    mapH = BitArray(undef,nr)
    mapH .= true
    mapVown = BitArray(undef,nr)
    mapVint = BitArray(undef,nr)
    
    for i = Int32(0):dH
        mapH .&= (hBoundInds[:,1].+i .<= hBoundInds[:,2])
        mapVown.=mapH
        for j = Int32(0):dVown
            mapVown .&= (vownBoundInds[:,1].+j .<= vownBoundInds[:,2])
            mapVint.=mapVown
            for k = Int32(0):dVint
                mapVint .&= (vintBoundInds[:,1].+k .<= vintBoundInds[:,2])
                append!(c,inds[mapVint])
                append!(r,indicesToIndex_Vector(vintBoundInds[mapVint,1].+k,vownBoundInds[mapVint,1].+j,hBoundInds[mapVint,1].+i))
            end
        end
    end
    sortIdx = sortperm(c)
    return r[sortIdx],c[sortIdx]
end
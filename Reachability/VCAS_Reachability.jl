using Printf
using HDF5
using Mmap
using GridInterpolations

include("NNet_Calculations.jl")
include("VCAS_Plotting.jl")
include("Reach_Functions.jl")
include("Reach_Constants.jl")


function getNetworks(folder="networks", ver=1, hu=25, epochs=1000)
    #= Load all of the neural networks
    Given information about neural network version, load the networks into a list
    
    NOTE: File format assumed here that may need to be changed
    Inputs: 
        folder (string): Folder where neural networks are stored
        ver (int, 2): version number of networks
        hu  (int, 40): Number of hidden units in each layer (written in filename)
        epochs (int, 1000): Epoch number of neural network
    
    Outputs:
        nets: List of NNet objects
    =#
    nets = []
    for pra in ACTIONS
        nnet_file_name = @sprintf("%s/VertCAS_newTransIntrSpeed_pra%02d_v%d_%dHU_%03d.nnet", folder, pra+1, ver, hu, epochs)
        nets = vcat(nets, NNet(nnet_file_name))
    end
    return nets
end

function getAdvParams(ra, delta1, delta2)
    cocAccelMin = -delta1 #-3-delta #-G/8 - delta
    cocAccelMax =  delta1 #3+delta #G/8 + delta
    weakAccelMin =  12.2-delta1 #-G/3 - delta
    weakAccelMax =  12.2+delta1 # G/4 + delta
    strongAccelMin =  13.4-delta1 #-G/3 - delta
    strongAccelMax =  13.4+delta1  # G/3 + delta
    
    params = Dict()
    params["aIntLo"] = -delta2
    params["aIntHi"] = delta2
    
    if ra == COC
        params["alo"] = cocAccelMin
        params["vlo"] = 100 #41.666
        params["ahi"] = cocAccelMax
        params["vhi"] = -100 #-41.666
    elseif ra == DNC
        params["alo"] = -weakAccelMax
        params["vlo"] = -100 #-41.666
        params["ahi"] = -weakAccelMin
        params["vhi"] = 0
    elseif ra == DND
        params["alo"] = weakAccelMin
        params["vlo"] = 0
        params["ahi"] = weakAccelMax
        params["vhi"] = 100 #41.666
    elseif ra == DES1500
        params["alo"] = -weakAccelMax
        params["vlo"] = -100 #-41.666 #100
        params["ahi"] = -weakAccelMin
        params["vhi"] = -25
    elseif ra == CL1500
        params["alo"] = weakAccelMin
        params["vlo"] = 25
        params["ahi"] = weakAccelMax
        params["vhi"] = 100 #41.666
    elseif ra == SDES1500
        params["alo"] = -strongAccelMax
        params["vlo"] = -100 #-41.666
        params["ahi"] = -strongAccelMin
        params["vhi"] = -25
    elseif ra == SCL1500
        params["alo"] = strongAccelMin
        params["vlo"] = 25
        params["ahi"] = strongAccelMax
        params["vhi"] = 100 #41.666
    elseif ra == SDES2500
        params["alo"] = -strongAccelMax
        params["vlo"] = -100 #-41.666
        params["ahi"] = -strongAccelMin
        params["vhi"] = -41.6666
    elseif ra == SCL2500
        params["alo"] = strongAccelMin
        params["vlo"] = 41.6666
        params["ahi"] = strongAccelMax
        params["vhi"] = 100 #41.666
    else
        println("UH OH")
    end
    return params
end

function getActionName(ra)
    #= Convert action index to a string name
    Inputs:
        ra (int): Advisory index
    Outputs:
        (string) name of advisory
    =#
    if ra == COC
        return "COC"
    elseif ra == DNC
        return "DNC"
    elseif ra == DND
        return "DND"
    elseif ra == DES1500
        return "DES1500"
    elseif ra == CL1500
        return "CL1500"
    elseif ra == SDES1500
        return "SDES1500"
    elseif ra == SCL1500
        return "SCL1500"
    elseif ra == SDES2500
        return "SDES2500"
    elseif ra == SCL2500
        return "SCL2500"
    else
        return "UNKNOWN RA"
    end
end

function indicesToIndex_Vector(vintInds, vownInds, hInds)
    #= Turn list of vint, vown, and h indices into a vector of cell indexes
    Get the cell index for each set of vint, vown, and h indices
    Inputs:
        vintInds (list of int): List of vint index values
        vownInds (list of int): List of vown index values
        hInds (list of int): List of h index values
    Outputs:
        List of cell indexes
    =#
    
    if all((vintInds.>=1) .& (vownInds.>=1) .& (hInds.>=1) .& (vintInds.<=NUMVINT) .& (vownInds.<=NUMVOWN) .& (hInds.<=NUMH))
        vownInds[vownInds.<Int32(1)] .+= NUMVOWN
        vintInds[vintInds.<Int32(1)] .+= NUMVINT
        return vintInds .+ NUMVINT.*(vownInds.-Int32(1)) .+ NUMVINT.*NUMVOWN.*(hInds.-Int32(1))
    end
    error("Bad Index")
end

function pointToIndex(h, vown, vint)
    #= Determine the cell that includes the given point. Truncate values if needed
    Inputs:
        h (float): h-value of point
        vown (float): vown-value of point
        vint (float): vint-value of point
    Outputs:
        Cell index
    =#
    if vint > VINTS[end-1]
        vint = VINTS[end-1]
    end
    if vown > VOWNS[end-1]
        vown = VOWNS[end-1]
    end
    if h > HS[end-1]
        h = HS[end-1]
    end
    hInd = findall(HS .> h)[1]-1
    vownInd = findall(VOWNS .> vown)[1]-1
    vintInd = findall(VINTS .> vint)[1]-1
    
    ind = indicesToIndex(vintInd,vownInd,hInd)
    return ind
end

function pointToSet(h, vown, vint)
    #= Given a point, return a one-hot vector of the cell indexes occupied by the point
    Inputs:
        h: h-value of point
        vown: vown-value of point
        vint: vint-value of point
    Outputs:
        One-hot vector of occupied cell
    =#
    ind = pointToIndex(h, vown, vint)
    set = zeros(Int32,NUMREGIONS)
    set[ind] = Int32(1)
    return set
end

function indicesToIndex(vintInd, vownInd, hInd)
    #= Turn vint, vown, and h indices into a cell index
    Inputs:
        pInd (int): Psi index value
        yInd (int): Y index value
        xInd (int): X index value
    Outputs:
        Cell index
    =#
    if vintInd <= NUMVINT &&  vownInd <= NUMVOWN && hInd <= NUMH && vintInd >= 1 && vownInd>=1 && hInd>=1
        return vintInd + NUMVINT*(vownInd-1) + NUMVINT*NUMVOWN*(hInd-1)
    end
    error("Bad Index")
end

function indexToIndices(ind)
    #= Turn a cell index into vint, vown, and h indices
    Inputs:
        Cell index (int)
    Outputs:
        vintInd (int): vint index value
        vownInd (int): vown index value
        hInd (int): h index value
    =#
    vintInd = mod(ind-1,NUMVINT)+1
    ind = div(ind-vintInd,NUMVINT)
    vownInd= mod(ind,NUMVOWN)+1
    hInd = div(ind-vownInd+1,NUMVOWN)+1
    return (vintInd,vownInd,hInd)
end

function indexToBounds(ind)
    #= Turn a cell index into the hyper-rectangular bounds of the cell
    Inputs:
        Cell index (int)
    Outputs:
        vint Lower (float)
        vint Upper (float)
        vown Lower (float)
        vown Upper (float)
        h Lower (float)
        h Upper (float)
    =#
    vintInd, vownInd, hInd = indexToIndices(ind)
    return (VINTS[vintInd],VINTS[vintInd+1],VOWNS[vownInd],VOWNS[vownInd+1],HS[hInd],HS[hInd+1])
end

function getAllBounds()
    #= Turn a cell index into vint, vown, and h indices
    Inputs:
        Cell index (int)
    Outputs:
        vintInd (int): vint index value
        vownInd (int): vown index value
        hInd (int): h index value
    =#
    hMin = [HS[k]    for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    hMax = [HS[k+1]  for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    vownMin = [VOWNS[j]   for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    vownMax = [VOWNS[j+1] for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    vintMin = [VINTS[i]   for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    vintMax = [VINTS[i+1] for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    return hcat(hMin,hMax), hcat(vownMin,vownMax), hcat(vintMin,vintMax)
end

function getAccelBounds(ra, delta1, delta2)
    #= Create arrays for the turning min/max values for each region
    Each region has the same turn rate bounds in this case
    Inputs:
        ra (int): resolution advisory index
        delta: (float): parameter to relax turn rate bounds
    Outputs:
        a0min (float array): minimum acceleration of ownship for each region
        a0max (float array): maximum acceleration of ownship for each region
        a1min (float array): minimum acceleration of intruder for each region
        a1max (float array): maximum acceleration of intruder for each region
    =#
    params = getAdvParams(ra, delta1, delta2)
    a0min = params["alo"]
    a0max = params["ahi"]
    a1min = params["aIntLo"]
    a1max = params["aIntHi"]
    v0min = params["vlo"]
    v0max = params["vhi"]
    return a0min,a0max,a1min,a1max,v0min,v0max
end


######################################################
### Functions for different initial reachable sets ###
######################################################

function getInitialSet(;pd=0)
    #= Create an initial reachable set. The keys are the sequences of most recent 
       advisories. Initially, the only key is a sequence of COC, since the system
       has not taken effect yet for the initial reachable set. For this function,
       all cells are reachable initially.
    Inputs:
        pd (int): Seconds of pilot delay. This determines the number of most recent 
            advisories that we need to track.
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    reachSet = Dict()
    initialKey = "0"*string(COC)
    for i = 2:pd
        initialKey *= string(COC)
    end
    reachSet[initialKey] = BitArray(undef,NUMREGIONS)
    reachSet[initialKey].= true
    return reachSet
end

function getInitialSet_point(pd=0; h=0, vown=0, vint=0)
    #= Same as getInitialSet(), except only the cell closest to the given point is
       included in the initial reachable set.
    Inputs:
        pd (int): Seconds of pilot delay. This determines the number of most recent 
            advisories that we need to track.
        h (float): h-value of point
        vown (float): vown-value of point
        vint (float): vint-value of point
    Outupts:
         reachSet (Dictionary of BitArrays): Return a dictionary that represents the 
            reachable set for different sequences of advisories. The dictionary entries
            are a BitArray of length NUMREGIONS, so each True in the array represents
            a cell that is reachable given that sequence of most recent advisories
    =#
    reachSet = Dict()
    initialKey = "0"*string(COC)
    for i = 1:pd
        initialKey *= string(COC)
    end
    index = pointToIndex(h, vown, vint)
    reachSet[initialKey] = BitArray(undef,NUMREGIONS)
    reachSet[initialKey].= false
    reachSet[initialKey][index] = true
    return reachSet
end

function getNmacCells()
    #= Compute cells that are NMACs (vertical separation less than 100 feet)
    Inputs: 
        None
    Outputs:
        nmacCells (BitArray): True for cells that are in the NMAC region of the state space
    =#
    hMin = [HS[k]    for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    hMax = [HS[k+1]  for i=1:NUMVINT, j=1:NUMVOWN, k=1:NUMH][:]
    return (hMin.<100) .& (hMax.>-100)
end
    

###############################################
### Helper functions for the reachable sets ###
###############################################

function saveSets(filename, sets)
    #= Write the sets dictionary to an h5py file
       The first key is the time
       The second key is the advisory sequence
    Inputs: 
        filename (string): Name of file to write
        sets (Dictionary): Reachable set dictionary
    =#
    h5open(filename,"w") do file
        for (k,v) in sets
            for (k2,v2) in v
                write(file,string(k)*"/"*k2,convert(Array{UInt8},v2))
            end
        end
    end
end

function loadSets(filename,maxT=40)
    #= Read the sets from an h5py to a dictionary
       The first key is the time
       The second key is the advisory sequence
    Inputs: 
        filename (string): Name of file to read
    Outputs:
        sets: 
    =#
    sets = Dict()
    h5open(filename,"r") do file
        for t in names(file)
            if parse(Int,t)<=maxT
                sets[parse(Int,t)] = Dict()
                for key in names(file[t])
                    sets[parse(Int,t)][key] = convert(BitArray,read(file[t][key]))
                end
            end
        end
    end
    return sets
end


function pruneSets!(sets,regKeep)
    #= Removes cells from the reachable set that are not in regKeep
    Inputs: 
        sets (Dictionary): Reachable set. The key is the sequence of previous advisories
        regKeep (BitArray): Array that is true if we should keep that cell in the set
    Outputs:
        None (The sets input is modified)
    =#
    for (k,v) = sets
        v[.!regKeep].=false
    end
end

function copySet(set)
    #= Deep copy a dictionary of BitArrays
    Inputs:
        set (Dictionary of BitArrays): Reachable set
    Outputs:
        out: Copy of reachable set
    =#
    out = Dict()
    for (k,v) = set
        out[k] = copy(v)
    end
    return out
end

function isConverged(sets,t;verbose=false)
    #= Check if a set of cells has converged at time t. Convergence
        is checked by seeing if any change has been made at the previous step.
        Since t counts down with each step, check if the reachable set at t+1
        is identical to the reachable set at time t.
    Inputs:
        Sets (Dictionary): Reachable sets dictionary. The first key is tau, 
            second key is most recent advisory sequence
        t (int): Tau value to check for convergence
    Outputs:
        allTrue (bool): True if the reachable set has not changed since the previous step
    
    =#
    allTrue = true
    sumDiff = 0
    if !(t+1 in keys(sets))
        return false
    end
    for key in keys(sets[t])
        if key in keys(sets[t+1])
            allTrue &= sets[t][key]==sets[t+1][key]
            if verbose
                sumDiff += sum(sets[t][key].!=sets[t+1][key])
            elseif !allTrue
                return false
            end
        else
            allTrue = false
        end
    end
    if verbose
        @printf("Number of different cells: %d\n",sumDiff)
    end
    return allTrue
end

function isNmac(sets,t,nmacCells)
    #= Check if any cells are reachable that are NMACS
       (near midair collisions)
    Inputs:
        Sets (Dictionary): Reachable sets dictionary. The first key is tau, 
            second key is most recent advisory sequence
        t (int): Tau value to check
        nmacCells (BitArray): True for cells that are NMACs
    Outputs:
        (bool): True if any cells in set are NMAC cells
    =#
    for (k,v) = sets[t]
        if any(v .& nmacCells)
            return true
        end
    end
    return false
end


##############################################################################
### Functions to create a memory-mapped array for over-approximated policy ###
##############################################################################

function sampleToMmap(nets)
    #= Evaluate network by sampling points at the center of each cell
    Inputs:
        nets (List of NNet): Neural networks to evaluate
    Outputs:
        acts (BitArray): True where the network could give an advisory in a cell
    =#
    acts = BitArray(undef,NUMTAU,NUMACTION,NUMACTION,NUMREGIONS)
    
    # Check center of cell. Change this to check a different point within each cell
    dh = 0.5; dvown = 0.5; dvint = 0.5; dtau = 0.5
    
    # Making Mesh
    println("Making Mesh")
    vintCenter = (1 - dvint) * VINTS[2:end] + dvint * VINTS[1:end-1]
    vownCenter = (1 - dvown) * VOWNS[2:end] + dvown * VOWNS[1:end-1]
    hCenter = (1 - dh) * HS[2:end] + dh * HS[1:end-1]
    hMesh = [h for vint=vintCenter, vown=vownCenter, h=hCenter]
    vownMesh = [vown for vint=vintCenter, vown=vownCenter, h=hCenter]
    vintMesh = [vint for vint=vintCenter, vown=vownCenter, h=hCenter]
    
    # Evaluate each network through batches
    # The output associated with the highest score is the advisory given by the network
    for (tauInd,tau) = enumerate(TAUS)
        tauMesh  = ones(NUMREGIONS)*tau 
        netIn = permutedims(hcat(hMesh[:],vownMesh[:],vintMesh[:],tauMesh),[2,1])
        for (praInd,pra) = enumerate(ACTIONS)
            @printf("Tau: %d, pRA: %d\n",tau,pra)
            net = nets[praInd]
            netOut = evaluate_network_multiple(net,netIn)
            _, mxindx = findmax(netOut,dims=1)
            for i = 1:length(mxindx)
                acts[tauInd,praInd,mxindx[i][1],i] = true
            end
        end
    end
    return acts
end

function reluValToMmap(folder,ver=5;praInds=1:NUMACTION)
    #= Read the output from ReluVal to create the network action BitArray
    Inputs:
        ver (int): Version of neural networks
        praInds (Array of int): previous advisories to read
        tauInds (Array of int): tau indices to read
    Outputs:
        acts (BitArray): True where the network could give an advisory in a cell
    =#
    acts = BitArray(undef,NUMTAU,NUMACTION,NUMACTION,NUMREGIONS)
    acts.=false
    
    for praInd = praInds
        @printf("pRA: %s\n",getActionName(ACTIONS[praInd]))
        file = @sprintf("%s/VC_newTransIntrSpeed_v%d_%02d.txt",folder,ver,praInd)
        raInd=1
        tauInd=0
        f=open(file)
        # SKIP FIRST LINE
        line=readline(f)
        @printf("\tRA: %s\n",getActionName(ACTIONS[raInd]))
        flush(stdout)
        while !eof(f)
            line=readline(f)
            if length(line)>0 && line[1]!='t' # Skip empty lines and time statements
                if line[1]=='R' # If the first letter is an 'R', advance to next ra index
                    tauInd+=1
                    if tauInd > NUMTAU
                        tauInd -= NUMTAU
                        raInd += 1
                        @printf("\tRA: %s\n",getActionName(ACTIONS[raInd]))
                        flush(stdout)
                    end
                    
                else
                    reg = parse(Int32,split(line,",")[1])+1 # +1 because Julia is 1-indexed
                    acts[tauInd,praInd,raInd,reg] = true; # Mark the tau/pra/ra/cell as true
                end
            end
        end
        close(f)
    end
    return acts
end

function writeNetworkActionsMmap(;folder="/raid/kjulian3/VertCAS/Reach/NetworkApprox",nnetFolder="networks",reluvalFolder="../ReluVal_intrSpeed",ver=1, ver2 = 1,hu=25, epochs=1000,useReluVal=false)
    #= Write advisories given by network to a memory-mapped array. Advisories can
       be determined by sampling a point within each cell or by using ReluVal. If 
       using ReluVal, then ReluVal must be run first.
    Inputs:
        folder (string): Folder to write memory-mapped file
        ver (int): Version of neural network
        hu (int): Hidden units in each layer of network
        epochs (int): Number of epochs used for training
        useReluVal (bool): True if using ReluVal, False for sampling
    Outputs:
        None (results are written to memory-mapped files)
    =#
    
    token="samp"
    if useReluVal
        token="reluVal"
    end
    
    filename=@sprintf("%s/netActions_mm_v%d.%d_%s_rect.bin",folder,ver,ver2,token)
    netActions=nothing
    if useReluVal
        @time netActions = reluValToMmap(reluvalFolder,ver)
    else
        nets = getNetworks(nnetFolder,ver,hu,epochs)
        @time netActions = sampleToMmap(nets);
    end
    
    # Write to memory mapped file. First write the dimensions of the array, then write the full array
    # Dimension of netActions are tau, previous advisory, next advisory, and cell index
    # The cells are 3D and have dimensions h, vown, and vint
    s = open(filename,"w+")
    write(s,Int64(size(netActions,1)))
    write(s,Int64(size(netActions,2)))
    write(s,Int64(size(netActions,3)))
    write(s,Int64(size(netActions,4)))
    write(s,netActions)
    close(s)
end

function getTables(tableFolder, ver)
    tables = Dict()
    for pra in ACTIONS
        tableFileName = @sprintf("%s/VertCAS_noResp_newTrans_TrainingData_v%d_%02d.h5", tableFolder, ver, pra+1)
        table = h5open(tableFileName, "r") do file
            read(file, "y")
        end;
        tables[pra] = table
    end
    return tables
end

function nearestNeighbor(grid, table, point)
    inds, weights = interpolants(grid, point)
    return argmax(table[:,inds[argmax(weights)]])
end 

function sampleTableToMmap(tables)
    #= Evaluate network by sampling points at the center of each cell
    Inputs:
        nets (List of NNet): Neural networks to evaluate
    Outputs:
        acts (BitArray): True where the network could give an advisory in a cell
    =#
    
    vels = vcat(LinRange(-100,-60,5),LinRange(-50,-35,4),LinRange(-30,30,21),LinRange(35,50,4),LinRange(60,100,5))
    relhs   = vcat(LinRange(-8000,-4000,5),LinRange(-3000,-1250,8),LinRange(-1000,-800,3),LinRange(-700,-150,12),LinRange(-100,100,9), LinRange(150,700,12),LinRange(800,1000,3),LinRange(1250,3000,8),LinRange(4000,8000,5))
    dh0s = vels
    dh1s = vels
    taus  = LinRange(0,40,41)
    grid = RectangleGrid(relhs,dh0s,dh1s,taus)
    
    acts = BitArray(undef,NUMTAU,NUMACTION,NUMACTION,NUMREGIONS)
    
    # Check center of cell. Change this to check a different point within each cell
    dhs = [0.0,1.0]
    dvowns = [0.0,1.0]
    dvints = [0.0,1.0]
    
    pointIndex = 1
    for dh in dhs
        for dvown in dvowns
            for dvint in dvints
                @printf("Point %d of %d\n",pointIndex,length(dhs)*length(dvowns)*length(dvints))
                vintCenter = (1 - dvint) * VINTS[2:end] + dvint * VINTS[1:end-1]
                vownCenter = (1 - dvown) * VOWNS[2:end] + dvown * VOWNS[1:end-1]
                hCenter = (1 - dh) * HS[2:end] + dh * HS[1:end-1]
                reg = 1
                percTarget = 1
                for h in hCenter
                    for vown in vownCenter
                        for vint in vintCenter
                            for (tauInd, tau) = enumerate(TAUS)
                                for (praInd, pra) = enumerate(ACTIONS)
                                    acts[tauInd, praInd, nearestNeighbor(grid, tables[pra],[h,vown,vint,tau]),reg] = true
                                end
                            end
                            reg+=1
                            if reg/NUMREGIONS >= percTarget/100.0
                                @printf("%d%% Complete\n",percTarget)
                                percTarget+=1
                            end
                        end
                    end
                end
                pointIndex+=1
            end
        end
    end
    return acts
end 
    
function writeTableActionsMmap(;folder="/raid/kjulian3/VertCAS/Reach/TableApprox",tableFolder="/raid/kjulian3/VertCAS/TrainingData",ver=5,ver2=4)
    #= Write advisories given by network to a memory-mapped array. Advisories can
       be determined by sampling a point within each cell or by using ReluVal. If 
       using ReluVal, then ReluVal must be run first.
    Inputs:
        folder (string): Folder to write memory-mapped file
        ver (int): Version of the table
        ver2 (int): Version of the state space discretization
    Outputs:
        None (results are written to memory-mapped files)
    =#

    filename=@sprintf("%s/tableActions_mm_v%d.%d_rect.bin",folder,ver,ver2)
    tables = getTables(tableFolder,ver)
    @time tableActions = sampleTableToMmap(tables);
    
    # Write to memory mapped file. First write the dimensions of the array, then write the full array
    # Dimension of tableActions are tau, previous advisory, next advisory, and cell index
    # The cells are 3D and have dimensions h, vown, and vint
    s = open(filename,"w+")
    write(s,Int64(size(tableActions,1)))
    write(s,Int64(size(tableActions,2)))
    write(s,Int64(size(tableActions,3)))
    write(s,Int64(size(tableActions,4)))
    write(s,tableActions)
    close(s)
end

function readNetworkActionsMmap(;folder="/raid/kjulian3/VertCAS/Reach/NetworkApprox",ver=5, ver2=4, token="samp")
    #= Read advisories given by network from a memory-mapped array.
    Inputs:
        folder (string): Folder to write memory-mapped file
        ver (int): Version of neural network
        token (string): Indentifying string (samp or reluVal)
    Outputs:
        netActions: Memory-mapped 4-D BitArray
    =#
    filename=@sprintf("%s/netActions_mm_v%d.%d_%s_rect.bin",folder,ver,ver2,token)
    s = open(filename)
    m1 = read(s,Int64)
    m2 = read(s,Int64)
    m3 = read(s,Int64)
    m4 = read(s,Int64)
    return Mmap.mmap(s,BitArray,(m1,m2,m3,m4))
end

function readTableActionsMmap(;folder="/raid/kjulian3/VertCAS/Reach/TableApprox",ver=5, ver2=4, token="samp")
    #= Read advisories given by network from a memory-mapped array.
    Inputs:
        folder (string): Folder to write memory-mapped file
        ver (int): Version of neural network
        token (string): Indentifying string (samp or reluVal)
    Outputs:
        netActions: Memory-mapped 4-D BitArray
    =#
    filename=@sprintf("%s/tableActions_mm_v%d.%d_rect.bin",folder,ver,ver2)
    s = open(filename)
    m1 = read(s,Int64)
    m2 = read(s,Int64)
    m3 = read(s,Int64)
    m4 = read(s,Int64)
    return Mmap.mmap(s,BitArray,(m1,m2,m3,m4))
end

###########################################################################
### Generate memory-mapped array of cell transitions as a sparse array  ###
###########################################################################

function writeReachDynamicsMmap(;folder="/raid/kjulian3/VertCAS/Reach/ReachDynamics",ras=ACTIONS,ver=1,delta1=0.0,delta2=0.0)
    #= Write reachableDynamics for each advisory to a memory-mapped array
    Inputs:
        folder (string): Folder to write memory-mapped file
        ras (list of int): List of advisories to consider
        delta (float): Parameter to relax the turn rate bounds
    Outputs:
        None (results are written to memory-mapped files)
    =#
    println("Writing reachable state dynamics to memory-mapped file")
    @printf("delta1 = %.3f, delta2 = %.3f\n",delta1,delta2)
    for ra = ras
        @printf("RA: %s\n",getActionName(ra))
        r,c = computeReachDynamics(delta1,delta2,ra)
        idx_ptr = [UInt32(1)]
        append!(idx_ptr,findall(c[1:end-1].!=c[2:end]).+UInt32(1))
        append!(idx_ptr,[UInt32(length(c)+1)])
        next_idx = UInt32.(r)

        s = open(@sprintf("%s/reachDynamics_mm_v%d_delta1%.2f_delta2%.2f_ra%d_idxPtr.bin",folder,ver,delta1,delta2,ra),"w+")
        write(s, UInt32(size(idx_ptr,1)))
        write(s, idx_ptr)
        close(s)

        s = open(@sprintf("%s/reachDynamics_mm_v%d_delta1%.2f_delta2%.2f_ra%d_nextIdx.bin",folder,ver,delta1,delta2,ra),"w+")
        write(s, UInt32(size(next_idx,1)))
        write(s, next_idx)
        close(s)
        r=nothing; c=nothing; idx_ptr=nothing; next_idx=nothing
    end
end

function readReachDynamicsMmap(;folder="/raid/kjulian3/VertCAS/Reach/ReachDynamics",ver=1,delta1=0.0,delta2=0.0,ras=ACTIONS)
    #= Read reachableDynamics for each advisory as a memory-mapped array
    Inputs:
        folder (string): Folder to read memory-mapped file
        ver (int): Version of dynamics
        delta (float): Parameter to relax the turn rate bounds
        ras (list of int): List of advisories to load
    Outputs:
        reachDynamics: List of memory-mapped vectors
    =#
    reachDynamics = Array{Vector{UInt32},1}[]
    
    for ra = ras
        fn = @sprintf("%s/reachDynamics_mm_v%d_delta1%.2f_delta2%.2f_ra%d_idxPtr.bin",folder,ver,delta1,delta2,ra)
        s = open(fn)
        m = read(s,UInt32)
        idxPtr = Mmap.mmap(s,Vector{UInt32},m)

        fn = @sprintf("%s/reachDynamics_mm_v%d_delta1%.2f_delta2%.2f_ra%d_nextIdx.bin",folder,ver,delta1,delta2,ra)
        s2 = open(fn)
        m = read(s2,UInt32)
        nextIdx = Mmap.mmap(s2,Array{UInt32,1},(m,));
    
        push!(reachDynamics,[idxPtr,nextIdx])
    end
    return reachDynamics
    
end

##########################################################
### Functions for computing the reachable set dynamics ###
##########################################################

function computeReachDynamics(delta1,delta2,ra)
    #= Compute the reachable set dynamics
    Inputs:
        delta (float): Paramter to relax turn rate bounds
        ra (int): Advisory
    Outputs:
        reachDynamics (list): Data structure representing which cells
            can be reached at the next time step from the current cell
    =#
    
    # Get the initial cell and turn rate bounds
    cellBounds = getAllBounds()
    accelBounds = getAccelBounds(ra,delta1,delta2)
    
    # Compute the region of the state space reachable at the next time step
    getNextSetBounds!(cellBounds, accelBounds)
    
    # Clear some memory
    accelBounds = nothing; GC.gc()
    
    # Compute the cells that overlap with the reachable region
    cellBoundInds = getIndexBounds(cellBounds)
    
    # Clear some memory
    cellBounds=nothing; GC.gc()
    
    # Return reachable dynamics array represented with vectors of pointers and indices
    return getReachDynamics(cellBoundInds)
end

############################################################################
### Function to compute the set of reachable cells at the next time step ###
############################################################################

function isReversal(pra,ra)
    if pra>DND && ra>COC && mod(pra,2)!=mod(ra,2)
        return true
    end
    return false
end

function getNextSetMmap(set,networkActions,reachDynamics,tau;pd=0,preventDoubleReversals=false, alwaysTrueRegions=nothing)
    #= Compute the next reachable set
    Inputs:
        set (Dictionary of BitArrays): Current reachable set
        networkActions (BitArray): Memory-mapped array that is True where the network 
            could give the advisory in that cell.
        reachDynamics (List): Memory-mapped data structure used to determine which cells are
            reachable at the next time step give a current cell and acting advisory
        tau: Tau value, used to determine which network is being used 
        alwaysTrueRegions (BitArray): True for a cell that should always be reachable, such
            as cells that are far away horizontally. An intruder could always appear on the horizon
    Outputs:
        nextSet (Dictionary of BitArrays): Next reachable set
    =#
    
    # Determine the index of tau so we use the correct neural network
    if tau<TAUS[1]
        tau=TAUS[1]
    end
    tauInd = findall(tau.>=TAUS)[end]
    
    nextSet = Dict()
    kl=1 #length of key
    for (key,currentCells) in set
        kl = length(key)
        
        delayAdvs = [parse(Int32,key[i]) for i=2:length(key)]   # Advisory the pilot could be following, delayed
        pRA       = parse(Int32,key[end]) # Advisory just given by the network at the previous time step
        if pd==0
            delayAdvs = []
        end
    
        # Loop through each possible next advisory
        for ra=ACTIONS 
            
            revYet = parse(Int32, key[1])
            realDelayAdvs = [ra]
            if length(delayAdvs)>0
                realDelayAdvs = vcat(delayAdvs,ra)
            end
            
            raActing = ra
            if isReversal(pRA, ra)
                if revYet == 1 && preventDoubleReversals
                    raActing = pRA
                    realDelayAdvs[end] = raActing
                end
                revYet = 1
            elseif ra == COC
                revYet = 0
            end
            
            # Get the next key and determine in which cells the advisory could be given
            cellsGiven = currentCells .& networkActions[tauInd,pRA+1,ra+1,:]
            if any(cellsGiven)
                nextKey = string(revYet)*key[3:end]*string(raActing)
                
                # Initialize reachable array if key is not yet in the dictionary
                if !(nextKey in keys(nextSet))
                    nextSet[nextKey] = BitArray(undef,NUMREGIONS)
                    nextSet[nextKey].=false
                end

                # Use reachDynamics to determine add cells to the next reachable set
                cg = findall(cellsGiven)
                for actingRA in realDelayAdvs
                    idxPtrs = reachDynamics[actingRA+1][1]
                    nextIdxs = reachDynamics[actingRA+1][2]
                    ns = nextSet[nextKey]

                    for cell = cg
                        ns[nextIdxs[idxPtrs[cell]:(idxPtrs[cell+1]-1)]].=true
                    end
                end
            end
        end
    end
    
    # Add cells that are always reachable
    if alwaysTrueRegions != nothing
        key = ""
        for i=1:kl
            key*=string(COC)
        end
        if !(key in keys(nextSet))
            nextSet[key] = BitArray(undef,NUMREGIONS)
            nextSet[key].= false
        end
        nextSet[key][alwaysTrueRegions] .= true   
    end
    return nextSet
end

function runReachability(delta1,delta2,pd;ver=5,ver2=1,maxTime=40,minTime=-50,token="samp",initPoint=nothing, preventDoubleReversals=false,stopAtNMAC=false,useTable=false,sets=nothing)
    # Read in memory-mapped files 
    if useTable
        acts = readTableActionsMmap(folder="/raid/kjulian3/VertCAS/Reach/TableApprox",ver=ver)
    else
        acts = readNetworkActionsMmap(folder="/raid/kjulian3/VertCAS/Reach/NetworkApprox",ver=ver,token=token); 
    end
    reachDynamics = readReachDynamicsMmap(folder="/raid/kjulian3/VertCAS/Reach/ReachDynamics",ver=ver2,delta1=delta1,delta2=delta2);
    
    # Dictionary for storing our reachable sets
    if sets == nothing
        sets = Dict(); 
        if initPoint == nothing
            sets[maxTime] = getInitialSet(pd=pd)
        else
            sets[maxTime] = getInitialSet_point(pd,h=initPoint[0],vown=initPoint[1],vint=initPoint[2])
        end
    else
        maxTime = minimum(keys(sets))
    end
    
    alwaysTrueRegions = nothing
    nmacCells = getNmacCells();
    
    # Compute reachable sets iteratively. 
    # Stop if convergence, an NMAC, or minTime is reached.
    foundNMAC = false
    for t = maxTime-1:-1:minTime
        # Compute next reachable set
        @time sets[t] = getNextSetMmap(sets[t+1],acts,reachDynamics,t,pd=pd,alwaysTrueRegions=alwaysTrueRegions,preventDoubleReversals=preventDoubleReversals)
        flush(stdout)
        if mod(t,10)==0
            @printf("\nt = %d\n",t)
        end

        # Unsafe NMAC defined as reaching an NMAC once tau reaches 0 (t<=0)
        if (t<=0) && (isNmac(sets,t,nmacCells))
            println("We reached an NMAC!")
            foundNMAC = true
            flush(stdout)
            if stopAtNMAC
                break
            end
        end

        # Check if the previous reachable set is the same as the current reachable set
        # If so, reachability has converged to a steady state solution
        if (t<0) && isConverged(sets,t,verbose=true)
            println("Reached steady-state!")
            break
        end
    end
    if foundNMAC
        println("Found NMAC")
    else
        println("No NMACs Found")
    end
    return sets
end

function computeMargin(sets)
    allReach = BitArray(undef,NUMREGIONS)
    allReach .= false
    for (key,reachVec) in sets[0]
        allReach .|= reachVec
    end

    margin = 1000.0
    for reg in findall(allReach)
        (_,_,_,_,h0,h1) = indexToBounds(reg)
        newMargin = min(abs(h0),abs(h1))
        if newMargin<margin
            margin=newMargin
        end
    end
    return margin
end

function loadReachSets(delta1,delta2,pd,useTable=false,preventDoubleReversals=false,maxT=40)
    folder = "/raid/kjulian3/VertCAS/Reach/JGCD_Sets"
    token = "Reluval"
    if useTable
        token = "Table"
    end
    pdr=""
    if preventDoubleReversals
        pdr="PDR"
    end
    filename = @sprintf("%s/Sets_%s_pd%d%s_delta1%.3f_delta2%.3f.h5",folder,token,pd,pdr,delta1,delta2)
    return loadSets(filename, maxT)
end

function saveReachability(delta1,delta2,pd,preventDoubleReversals=false,sets=nothing,stopAtNMAC=true)
    folder = "/raid/kjulian3/VertCAS/Reach/JGCD_Sets"
    token = "reluVal"
    pdr=""
    if preventDoubleReversals
        pdr="PDR"
    end
    filename = @sprintf("Sets_Reluval_pd%d%s_delta1%.3f_delta2%.3f.h5",pd,pdr,delta1,delta2)
    println("\n"*filename)
    sets = runReachability(delta1,delta2,pd,token=token,stopAtNMAC=stopAtNMAC,useTable=false,preventDoubleReversals=preventDoubleReversals,sets=sets);
    saveSets(folder*"/"*filename, sets)
    
    #filename = @sprintf("Sets_Table_pd%d%s_delta1%.3f_delta2%.3f.h5",pd,pdr,delta1,delta2)
    #println(filename)
    #sets = runReachability(delta1,delta2,pd,token=token,stopAtNMAC=true,useTable=true,preventDoubleReversals=preventDoubleReversals);
    #saveSets(folder*"/"*filename, sets);
end
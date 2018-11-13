#Julia implementation of "load_network" function
export NNet, evaluate_network, evaluate_network_multiple, num_inputs
type NNet
    file::AbstractString
    weights::Array{Any,1}
    biases::Array{Any,1}
    symmetric::Int32
    numLayers::Int32
    inputSize::Int32
    outputSize::Int32
    maxLayerSize::Int32
    
    layerSizes::Array{Int32,1}
    mins::Array{Float64,1}
    maxes::Array{Float64,1}
    means::Array{Float64,1}
    ranges::Array{Float64,1}
    
    numNetworks::Int32
    numTau::Int32
    numPra::Int32
    tauCut::Array{Float64,1}
    praCut::Array{Float64,1}
    networkArray::Bool
    
    function NNet(file::AbstractString)
        this  = new()
        this.file = file
        f = open(this.file)
        
        ###############
        # Begin reading .nnet file
        ###############
        
        #Read header line
        line = readline(f)
        this.networkArray = false
        if line[1] == 'A'
            this.networkArray = true
        end
        
        if this.networkArray
            #Read number of neural networks
            line = readline(f)
            record = split(line,[',','\n'])
            this.numNetworks = parse(Int32,record[1])
            
            #Read number of previous RA cutpoints
            line = readline(f)
            record = split(line,[',','\n'])
            this.numPra = parse(Int32,record[1])
            
            #Read each previous RA cutpoint
            line = readline(f)
            record = split(line,[',','\n'])
            this.praCut = zeros(this.numPra)
            for i=1:(this.numPra)
                this.praCut[i] = parse(Float64,record[i])
            end
            
            #Read the number of tau cutpoints
            line = readline(f)
            record = split(line,[',','\n'])
            this.numTau = parse(Int32,record[1])
            
            #Read each tau cutpoint
            line = readline(f)
            record = split(line,[',','\n'])
            this.tauCut = zeros(this.numTau)
            for i=1:(this.numTau)
                this.tauCut[i] = parse(Float64,record[i])
            end
                
        else
            this.numNetworks = 1
        end
        
        #Read information about the neural network
        line = readline(f)
        record = split(line,[',','\n'])
        this.numLayers = parse(Int32,record[1])
        this.inputSize = parse(Int32,record[2])
        this.outputSize = parse(Int32,record[3])
        this.maxLayerSize=parse(Int32,record[4])
        
        line = readline(f)
        record = split(line,[',','\n'])
        this.layerSizes = zeros(this.numLayers+1)
        for i=1:(this.numLayers+1)
            this.layerSizes[i]=parse(Int32,record[i])
        end
        
        line = readline(f)
        record = split(line,[',','\n'])
        this.symmetric = parse(Int32,record[1])
        
        line = readline(f)
        record = split(line,[',','\n'])
        this.mins = zeros(this.inputSize)
        for i=1:(this.inputSize)
            this.mins[i]=parse(Float64,record[i])
        end
        
        line = readline(f)
        record = split(line,[',','\n'])
        this.maxes = zeros(this.inputSize)
        for i=1:(this.inputSize)
            this.maxes[i]=parse(Float64,record[i])
        end
        
        
        line = readline(f)
        record = split(line,[',','\n'])
        this.means = zeros(this.inputSize+1)
        for i=1:(this.inputSize+1)
            this.means[i]=parse(Float64,record[i])
        end
        
        line = readline(f)
        record = split(line,[',','\n'])
        this.ranges = zeros(this.inputSize+1)
        for i=1:(this.inputSize+1)
            this.ranges[i]=parse(Float64,record[i])
        end
        
        this.weights = Any[zeros(this.layerSizes[2],this.layerSizes[1])]
        this.biases  = Any[zeros(this.layerSizes[2])]
        for i=2:this.numLayers
            this.weights = [this.weights;Any[zeros(this.layerSizes[i+1],this.layerSizes[i])]]
            this.biases  = [this.biases;Any[zeros(this.layerSizes[i+1])]]
        end
        this.weights = Any[this.weights]
        this.biases  = Any[this.biases]
        for i=2:this.numNetworks
            w2 = Any[zeros(this.layerSizes[2],this.layerSizes[1])]
            b2 = Any[zeros(this.layerSizes[2])]
            for i=2:this.numLayers
                w2 = [w2;Any[zeros(this.layerSizes[i+1],this.layerSizes[i])]]
                b2 = [b2;Any[zeros(this.layerSizes[i+1])]]
            end
            this.weights = [this.weights;Any[w2]]
            this.biases  = [this.biases;Any[b2]]
        end
        
        nnetInd = 1
        layer=1
        i=1
        j=1
        line = readline(f)
        record = split(line,[',','\n'])
        while !eof(f)
            if layer>this.numLayers
                layer=1
                nnetInd=nnetInd+1
            end
            while i<=this.layerSizes[layer+1]
                while record[j]!=""
                    this.weights[nnetInd][layer][i,j] = parse(Float64,record[j])
                    j=j+1
                end
                j=1
                i=i+1
                line = readline(f)
                record = split(line,[',','\n'])
            end
            i=1
            while i<=this.layerSizes[layer+1]
                this.biases[nnetInd][layer][i] = parse(Float64,record[1])
                i=i+1
                line = readline(f)
                record = split(line,[',','\n'])
            end
            layer=layer+1
            i=1
            j=1
        end
        close(f)
        
        return this
    end
end

#Evaluates one set of inputs
function evaluate_network(nnet::NNet,input::Array{Float64,1})
    numLayers = nnet.numLayers
    inputSize = nnet.inputSize
    outputSize = nnet.outputSize
    symmetric = nnet.symmetric
    biases = nnet.biases
    weights = nnet.weights
    
    nnetInd = 1
    if nnet.numNetworks>1
        tau = input[6]
        pra = input[7]
        nnetInd += round(Int32,pra) - 2
        nnetInd *= nnet.numTau
        
        upper = nnet.numTau
        lower = 1
        while lower < upper
            middle = floor(Int32,(upper+lower)/2)
            if tau - nnet.tauCut[middle] > nnet.tauCut[middle+1]-tau
                lower = middle + 1
            else
                upper = middle
            end
        end
        nnetInd += lower
    end
    
    inputs = zeros(inputSize)
    for i = 1:inputSize
        if input[i]<nnet.mins[i]
            inputs[i] = (nnet.mins[i]-nnet.means[i])/nnet.ranges[i]
        elseif input[i] > nnet.maxes[i]
            inputs[i] = (nnet.maxes[i]-nnet.means[i])/nnet.ranges[i] 
        else
            inputs[i] = (input[i]-nnet.means[i])/nnet.ranges[i] 
        end
    end

    for layer = 1:numLayers-1
        temp = max.(*(weights[nnetInd][layer],inputs[1:nnet.layerSizes[layer]])+biases[nnetInd][layer],0)
        inputs = temp
    end
    outputs = *(weights[nnetInd][end],inputs[1:nnet.layerSizes[end-1]])+biases[nnetInd][end]
    for i=1:outputSize
        outputs[i] = outputs[i]*nnet.ranges[end]+nnet.means[end]
    end
    return outputs
end

#Evaluates multiple inputs at once. Each set of inputs should be a column in the input array
#Returns a column of output Q values for each input set
function evaluate_network_multiple(nnet::NNet,input::Array{Float64,2})
    numLayers = nnet.numLayers
    inputSize = nnet.inputSize
    outputSize = nnet.outputSize
    symmetric = nnet.symmetric
    biases = nnet.biases
    weights = nnet.weights
        
    _,numInputs = size(input)
    symmetryVec = zeros(numInputs)
    
    
    nnetInd = 1
    if nnet.numNetworks>1
        tau = input[6]
        pra = input[7]
        nnetInd += round(Int32,pra) - 2
        nnetInd *= nnet.numTau
        
        upper = nnet.numTau
        lower = 1
        while lower < upper
            middle = floor(Int32,(upper+lower)/2)
            if tau - nnet.tauCut[middle] > nnet.tauCut[middle+1]-tau
                lower = middle + 1
            else
                upper = middle
            end
        end
        nnetInd += lower
    end
    
    inputs = zeros(inputSize,numInputs)
    for i = 1:inputSize
        for j = 1:numInputs
            if input[i,j]<nnet.mins[i]
                inputs[i,j] = (nnet.mins[i]-nnet.means[i])/nnet.ranges[i]
            elseif input[i,j] > nnet.maxes[i]
                inputs[i,j] = (nnet.maxes[i]-nnet.means[i])/nnet.ranges[i] 
            else
                inputs[i,j] = (input[i,j]-nnet.means[i])/nnet.ranges[i] 
            end
        end
    end

    for layer = 1:numLayers-1
        inputs = max.(*(weights[nnetInd][layer],inputs[1:nnet.layerSizes[layer],:])+*(biases[nnetInd][layer],ones(1,numInputs)),0)
    end
    outputs = *(weights[nnetInd][end],inputs[1:nnet.layerSizes[end-1],:])+*(biases[nnetInd][end],ones(1,numInputs))
    for i=1:outputSize
        for j=1:numInputs
            outputs[i,j] = outputs[i,j]*nnet.ranges[end]+nnet.means[end]
        end
    end
    return outputs
end

function num_inputs(nnet::NNet)
    if nnet.numNetworks==1
        return nnet.inputSize
    else
        return nnet.inputSize+2
    end
end
function num_outputs(nnet::NNet)
    return nnet.outputSize
end
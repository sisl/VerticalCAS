export viz_policy


"""
Helper function to draw aircraft
"""
function getACString(theta,x,y,fill,draw,width=2.0)
    return "\\node[aircraft side,fill="*fill*",draw="*draw*", minimum width="*string(width)*"cm,rotate="*string(theta)*",scale = 0.35] at (axis cs:"*string(x)*", "*string(y)*") {};"
end


function get_belief(pstate::Vector{Float64}, grid::RectangleGrid,interp::Bool=false)
    belief = spzeros(nstates, 1)
    indices, weights = interpolants(grid, pstate)
    if !interp
        indices = indices[findmax(weights)[2]]
        weights = 1.0
    end
    for i = 1:length(indices)
        belief[indices[i]] = weights[i]
    end # for i
    return belief
end # function get_belief

function getRevActions(pra)
    if (pra>1) && (mod(pra,2)==0)
        return [3,5,7,9]
    elseif (pra>1)
        return [2,4,6,8]
    end
    return []
end

function prevRev(qvals,preventFlag,pra)
   if preventFlag
       for ra in getRevActions(pra)
            qvals[ra]-=1e8
        end
    end
    return qvals
end

function viz_policy(nnet_folder::AbstractString="",nnet_prefix::AbstractString="",table_folder::AbstractString="",table_prefix::AbstractString="",vnum::Int64=4, hu::Int64=45, epoch::Int64=200)
    
    # Make sure one of the folder paths is given
    if nnet_folder=="" && table_folder==""
        println("Need at least either a valid nnet or table folder path!")
        return
    end
    
    # Important variables
    currentPra = -1
    nnet = []
    Q = nothing
    grid = RectangleGrid(relhs,dh0s,dh1s,taus)
    
    
    
    COC = RGB(1.,1.,1.) # white
    SR = RGB(.0,.0,.5) # navy
    SL = RGB(.0,.600,.0) # green
    WR = RGB(.5,.5,.5) # grey
    WL = RGB(.7,.9,.0) # neon green
    ra_colors = [SL,WL,COC,WR,SR]
    bg_colors = [COC]
    
    # Colors
    ra_1 = RGB(1.,1.,1.) # white
    ra_2 = RGB(.0,1.0,1.0) # cyan
    ra_3 = RGB(144.0/255.0,238.0/255.0,144.0/255.0) # lightgreen
    ra_4 = RGB(30.0/255.0,144.0/255.0,1.0) # dodgerblue
    ra_5 = RGB(0.0,1.0,.0) # lime
    ra_6 = RGB(0.0,0.0,1.0) # blue
    ra_7 = RGB(34.0/255.0,139.0/255.0,34.0/255.0) # forestgreen
    ra_8 = RGB(0.0,0.0,128.0/255.0) # navy
    ra_9 = RGB(0.0,100.0/255.0,0.0) # darkgreen
    
    
    ## FOR PAPER
    ra_2 = WR
    ra_4 = SR #RGB(0.0,0.0,1.0)
    ra_6 = SR
    ra_8 = RGB(0.2, 0.0, 0.3)
    
    ra_3 = WL
    ra_5 = SL #RGB(0.0,1.0,0.0)
    ra_7 = SL
    ra_9 = RGB(.2,.300,.0)
    
    
    colors = [ra_1;ra_2;ra_3;ra_4;ra_5;ra_6;ra_7;ra_8;ra_9]
    bg_colors = [RGB(1.0,1.0,1.0)]
    
    # Create scatter plot classes for color key
    sc_string = "{"
    for i=1:9
        define_color("ra_$i",  colors[i])
        if i==1
            sc_string *= "ra_$i={mark=square, style={black, mark options={fill=ra_$i}, mark size=6}},"
        else
            sc_string *= "ra_$i={style={ra_$i, mark size=6}},"
        end
    end
    
    # Color key as a scatter plot
    sc_string=sc_string[1:end-1]*"}"
    xx = [-1.5,-1.5,-1.5, -1.5, -1.5, -1.5, -1.5, 0.4 ,.4,]
    yy = [1.65,1.15,0.65, 0.15, -0.35, -0.85, -1.35, 1.65, 1.15]
    zz = ["ra_1","ra_2","ra_3","ra_4","ra_5","ra_6","ra_7","ra_8","ra_9"]
    sc = string(sc_string)

    # Set up interactive display
    @manipulate for nbin = 100,
        savePlot = [false,true],
        xmin = 0.0,
        xmax = 40.0,
        ymin = -1000.0,
        ymax = 1000.0,
        dh0 = 0.0, 
        dh1 = 0.0, 
        pra = action_names,
        preventRev = [false,true]
        
        # Get previous RA index
        pra = findall(pra.==action_names)[1]
        if pra != currentPra
            if table_folder !=""
                Q = h5open(@sprintf("%s/%s_%02d.h5",table_folder,table_prefix,pra), "r") do file
                    read(file, "y")
                end
            end
            
            if nnet_folder != "" 
                if epoch>0
                    nnet = read_network(@sprintf("%s/%s_pra%02d_v%d_%dHU_%03d.nnet",nnet_folder,nnet_prefix,pra,vnum,hu,epoch))
                else
                    nnet = read_network(@sprintf("%s/%s_pra%02d_v%d_%dHU.nnet",nnet_folder,nnet_prefix,pra,vnum, hu))
                end
            end
            currentPra = pra
        end
        
        # Evaluate network at pixel points
        if nnet_folder != ""
            if num_inputs(nnet)==3
                inputsNet = hcat([[relh,dh0,tau] for tau=LinRange(xmin, xmax,nbin) for relh=LinRange(ymin,ymax,nbin)]...)
                q_nnet = evaluate_network_multiple(nnet,inputsNet)
            else        
                inputsNet = hcat([[relh,dh0,dh1,tau] for tau=LinRange(xmin, xmax,nbin) for relh=LinRange(ymin,ymax,nbin)]...)
                q_nnet = evaluate_network_multiple(nnet,inputsNet)
            end
        end
        
        # Q Table Heat Map
        function get_heat1(x::Float64, y::Float64)
            tau = x 
            relh = y
            bel = get_belief([relh,dh0,dh1,tau],grid,false)
            qvals = Q[:,bel.rowval[1]]
            qvals = prevRev(qvals,preventRev,pra)                                   
            return actions[findmax(qvals)[2]]
        end # function get_heat1
        
        ind = 1
        #Neural Net Heat Map
        function get_heat2(x::Float64, y::Float64)                        
            qvals = q_nnet[:,ind]
            qvals = prevRev(qvals,preventRev,pra)  
            ind +=1
            return actions[findmax(qvals)[2]]
        end # function get_heat2
        
        #Plot table or neural network policies if possible
        g = GroupPlot(2, 2, groupStyle = "horizontal sep=4cm, vertical sep=2.5cm")
        if table_folder!=""
            push!(g, Axis([
                Plots.Image(get_heat1, (xmin, xmax), (ymin, ymax), zmin = 1, zmax = 9,
                xbins = nbin, ybins = nbin, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false),
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Original Table Advisories"))
        end
        if nnet_folder!=""
            push!(g, Axis([
                Plots.Image(get_heat2, (xmin, xmax), (ymin, ymax), zmin = 1, zmax = 9, 
                            xbins = nbin, ybins = nbin, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false), 
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax,width="10cm", height="8cm", 
                   xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Neural Network Advisories"))
        end
             
        # Save policy to a tex file to be used in papers
        if savePlot
            PGFPlots.save("PolicyPlot.tex", g, include_preamble=true)
        else
                            
            # Create Color Key
            f = (x,y)->x # Dummy function for background white image
            push!(g, Axis([
                Plots.Image(f, (-2,2), (-2,2),colormap = ColorMaps.RGBArrayMap(bg_colors),colorbar=false),
                Plots.Scatter(xx, yy, zz, scatterClasses=sc),
                Plots.Node("RA 1: COC ",0.15,0.915,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 2: DNC ",0.15,0.790,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 3: DND",0.15,0.665,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 4: DES15000",0.15,0.540,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 5: CL1500 ",0.15,0.415,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 6: SDES1500",0.15,0.290,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 7: SCL1500",0.15,0.165,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 8:  SDES2500",0.63,0.915,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 9:  SCL2500",0.63,0.790,style="black,anchor=west", axis="axis description cs"),
                ],width="10cm",height="8cm", hideAxis =true, title="KEY"))
        end
    end
end 
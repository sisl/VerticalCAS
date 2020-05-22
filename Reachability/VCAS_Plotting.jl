using PGFPlots
using Reel
using Colors
using LaTeXStrings

# Setup plotting for animinations
# Set default axis and label styles
Reel.extension(m::MIME"image/svg+xml") = "svg"
Reel.set_output_type("gif") # may be necessary for use in IJulia


function plotReachable(sets,extent,time;zmax=50,timeOffset=0,cmCustom=true,desiredKey=nothing,vint=nothing,nbin=150,width="18cm",height="18cm",title=-1,xlabel=L"$\dot{h}_\text{own}\;(\text{kft/min})$",ylabel=L"$h\;(\text{kft})$",filename="Pics/VCAS_ReachablePlot.png",savePlot=false,saveName=saveName="Pics/VCAS_ReachablePlot.tex",isRow=false,isLast=false,keepReg=nothing)
    #= Plot reachable set
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        time (int): Time (tau) value to plot (Should be a key in sets)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom Jet colormap
        desiredKey (string): Previous advisory sequence to plot (should be a key of sets[time]).
            If nothing, the reachable sets of all previous advisory sequences are combined
        psi (float): Show reachable cells that contain a given psi value.
            If nothing, all reachable cells are used
        nbin (int): Number of bins for plotting resolution
        width (string): Width of plot
        height (string): height of plot
        title (string): Title of plot. If title==-1, use a default title
        xlabel (string): X-axis label
        ylabel (string): Y-axis label
        filename (string): Filename for saving png file (if saving)
        savePlot (bool): True if saving the plot as .tex and .png files
        saveName (string): Name for the .tex file
        isRow (bool): True if this plot is a part of a larger row of plots
        isLast (bool): True if this is the last plot in an animation
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        g (PGFPlots.Axis): Plot of reachable regions
    =#
    
    set = BitArray(undef,NUMREGIONS)
    if desiredKey != nothing
        set .= sets[time-timeOffset][desiredKey]
    else
        set .= false
        for (k,v) = sets[time-timeOffset]
            set .|= v
        end
    end
    
    if !cmCustom
        zmax=1
    end
    
    function getXY(x,y)
        y*=1000
        x*=1000.0/60.0
        
        
        if y>HS[end-1]
            y=HS[end-1]
        elseif y<HS[1]
            y=HS[1]
        end
        if x>VOWNS[end-1]
            x=VOWNS[end-1]
        elseif x<VOWNS[1]
            x=VOWNS[1]
        end
        
        xInd = findall(VOWNS.>x)[1]-1
        yInd = findall(HS.>y)[1]-1
        sumFound=0
        if vint !=nothing
            if vint>VINTS[end-1]
                vint=VINTS[end-1]
            elseif vint<VINTS[1]
                vint = VINTS[1]
            end
            vintInd = findall(VINTS.>vint)[1]-1
            ind = indicesToIndex(vintInd,xInd,yInd)
            if (keepReg != nothing) && (!keepReg[ind])
                return zmax
            end
            sumFound = set[ind]
        else
            inds = indicesToIndex_Vector(1:NUMVINT,xInd*ones(Int32,NUMVINT),yInd*ones(Int32,NUMVINT))
            if (keepReg != nothing) && (!keepReg[inds[1]])
                return zmax
            end
            sumFound = sum([set[i] for i in inds])
        end
        return sumFound
    end
    
    if !savePlot && !isRow
        filename=nothing
    end
    
    
    PGFPlots.resetPGFPlotsPreamble()
    PGFPlots.pushPGFPlotsPreamble("\\usepackage{amsmath}")
    if !isRow
        PGFPlots.pushPGFPlotsPreamble("\\pgfplotsset{every axis/.append style={title style={font=\\huge},label style={font=\\huge},tick label style={font=\\Large}}}")
    end
        
    xmin, xmax, ymin, ymax = extent
    if title==-1
        title=@sprintf("Time to loss of horizontal separation: %d s",max(0,time))
        if time<0
            title=@sprintf("Time since loss of horizontal separation: %d s",-time)
        elseif time==0
            title = "Lost horizontal separation" 
        end
        if isLast
            if time < 0
                title="Converged to steady state"
            else
                title="Reached NMAC"
            end
        end
    end
    if isRow
        title = latexstring("\\tau=",time,"\\text{s}")
        if isLast
            title="Converged"
        end
    end
    
    style=""
    if ylabel==""
        style = "yticklabels={,,}, scaled y ticks=false"
    end
    cm = ColorMaps.Named("Blues")
    if cmCustom
        viridis = [ "#ffffffff", "#440558ff", "#450a5cff", "#450e60ff", "#451465ff", "#461969ff",
                    "#461d6dff", "#462372ff", "#472775ff", "#472c7aff", "#46307cff", "#45337dff",
                    "#433880ff", "#423c81ff", "#404184ff", "#3f4686ff", "#3d4a88ff", "#3c4f8aff",
                    "#3b518bff", "#39558bff", "#37598cff", "#365c8cff", "#34608cff", "#33638dff",
                    "#31678dff", "#2f6b8dff", "#2d6e8eff", "#2c718eff", "#2b748eff", "#29788eff",
                    "#287c8eff", "#277f8eff", "#25848dff", "#24878dff", "#238b8dff", "#218f8dff",
                    "#21918dff", "#22958bff", "#23988aff", "#239b89ff", "#249f87ff", "#25a186ff",
                    "#25a584ff", "#26a883ff", "#27ab82ff", "#29ae80ff", "#2eb17dff", "#35b479ff",
                    "#3cb875ff", "#42bb72ff", "#49be6eff", "#4ec16bff", "#55c467ff", "#5cc863ff",
                    "#61c960ff", "#6bcc5aff", "#72ce55ff", "#7cd04fff", "#85d349ff", "#8dd544ff",
                    "#97d73eff", "#9ed93aff", "#a8db34ff", "#b0dd31ff", "#b8de30ff", "#c3df2eff",
                    "#cbe02dff", "#d6e22bff", "#e1e329ff", "#eae428ff", "#f5e626ff", "#fde725ff"]
        
        cm=ColorMaps.RGBArrayMap([convert(RGB{Float64},parse(Colorant,vir)) for vir in viridis], interpolation_levels=100)
    end
    colorbar=true
    colorbarStyle = "at={(1.03,0.07)},height=0.93*\\pgfkeysvalueof{/pgfplots/parent axis width},anchor=south west,ymin=1.0,ylabel=Number of reachable cells,ytick={1,5,10,15,20,25,30,35,40,45,50}"
    extraStyle = "colorbar/draw/.append code={\\begin{axis}[every colorbar,ytick={0},ylabel=,height=0.035*\\pgfkeysvalueof{/pgfplots/parent axis width},ymin=0,ymax=0,at={(1.03,0.0)},anchor=south west,colorbar=false]\\end{axis}}"
    
    if isRow
        if !isLast
            colorbar = false
            colorbarStyle=""
            extraStyle = ""
        else
            colorbarStyle = "at={(1.03,0.2)},height=2.27*\\pgfkeysvalueof{/pgfplots/parent axis width},anchor=south west,ymin=1.0,ylabel=Number of reachable cells,ytick={1,10,20,30,40,50}"
            
            colorbarStyle = "at={(1.03,0.13)},height=0.87*\\pgfkeysvalueof{/pgfplots/parent axis width},anchor=south west,ymin=1.0,ylabel=Number of reachable cells,ytick={1,10,20,30,40,50}"
        end
    end
    
    textPlots = Array{Plots.Node,1}()
    if isLast
        if isRow
            textPlots = [Plots.Node("Reachable",-0.5,0.8,style="black"),
                         Plots.Node("Unreachable", -0.5,0.3,style="black")]
        else
            textPlots = [Plots.Node("\\huge Reachable",0.6,1.200,style="black"),
                         Plots.Node("\\huge Unreachable", 0.6,0.3,style="black"),
                         Plots.Node("\\huge NMAC Region", -1.2,0.0,style="black")]
        end
    end
    
    g = Axis(vcat([
            Plots.Image(getXY,(xmin,xmax),(ymin,ymax),zmin=0,zmax=zmax,
                xbins=nbin,ybins=nbin,colormap=cm,colorbar=colorbar,colorbarStyle=colorbarStyle,filename=filename),
            Plots.Linear([-2.500,2.500],[-0.100,-.100],style="red!70!white,dashed,no marks,ultra thick"),
            Plots.Linear([-2.500,2.500],[0.100,.100],style="red!70!white,dashed,no marks, ultra thick")],
            textPlots),
        width=width,height=height,xlabel=xlabel,ylabel=ylabel,title=title,style=style*extraStyle)
    if savePlot
        save(saveName,g,include_preamble=false)
    else
        return g
    end
end


function plotReachableRow(sets,extent,times;nbin=100,zmax=50,timeOffset=0,cmCustom=true,savePlot=false,saveName="Pics/VCAS_ReachableRow.tex",keepReg=nothing)
    #= Plot a row of reachable plots, which are the reachable regions at different times
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        times (array of int): Times (taus) at which we want to plot (Should be keys in sets)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom viridis colormap
        nbin (int): Number of bins for plotting resolutiones
        saveName (string): Name for the .tex file
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        g (PGFPlots.GroupPlot): Row of reachable region plots
    =#
    PGFPlots.resetPGFPlotsPreamble()
    extraStyle = ",colorbar/draw/.append code={\\begin{axis}[every colorbar,ytick={0},ylabel=,height=.07*\\pgfkeysvalueof{/pgfplots/parent axis width},ymin=0,ymax=0,at={(1.03,0.0)},anchor=south west,colorbar=false]\\end{axis}}"
    g = GroupPlot(length(times),1, groupStyle = "horizontal sep=0.35cm", style="height=5.0cm, width=5.0cm"*extraStyle)
    width=nothing; height=nothing;
    xlabel=L"$\dot{h}_\text{own}\;(\text{kft/min})$";ylabel=L"$h\;(\text{kft})$"
    
    for time in times
        fn = nothing
        if savePlot
            fn = @sprintf("Pics/VCAS_ReachRow_%d.png",time)
        end
        push!(g,plotReachable(sets,extent,time,nbin=nbin,width=width,height=height,title=-1,xlabel=xlabel,ylabel=ylabel,filename=fn,zmax=zmax,timeOffset=timeOffset,cmCustom=cmCustom,isRow=true,keepReg=keepReg,isLast=time==times[end]))
        ylabel=""
    end
    if savePlot
        save(saveName,g,include_preamble=false)
        return
    else
        return g
    end
end

function plotReachableRow2(sets,extent,times;nbin=100,zmax=50,timeOffset=0,cmCustom=true,savePlot=false,saveName="Pics/VCAS_ReachableRow.tex",keepReg=nothing)
    #= Plot a row of reachable plots, which are the reachable regions at different times
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        times (array of int): Times (taus) at which we want to plot (Should be keys in sets)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom viridis colormap
        nbin (int): Number of bins for plotting resolutiones
        saveName (string): Name for the .tex file
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        g (PGFPlots.GroupPlot): Row of reachable region plots
    =#
    PGFPlots.resetPGFPlotsPreamble()
    extraStyle = ",colorbar/draw/.append code={\\begin{axis}[every colorbar,ytick={0},ylabel=,height=.12*\\pgfkeysvalueof{/pgfplots/parent axis width},ymin=0,ymax=0,at={(1.03,0.0)},anchor=south west,colorbar=false]\\end{axis}}"
    g = GroupPlot(round(Int,length(times)/2),2, groupStyle = "horizontal sep=0.35cm,vertical sep=1.5cm", style="height=4.8cm, width=4.8cm"*extraStyle)
    width=nothing; height=nothing;
    xlabel=L"$\dot{h}_\text{own}\;(\text{kft/min})$";ylabel=L"$h\;(\text{kft})$"
    
    for (i,time) in enumerate(times)
        fn = nothing
        if savePlot
            fn = @sprintf("Pics/VCAS_ReachRow_%d.png",time)
        end
        ylab=""
        xlab=""
        if mod(i,round(Int,length(times)/2))==1
            ylab=ylabel
        end
        if i>round(Int,length(times)/2)
            xlab=xlabel
        end
        push!(g,plotReachable(sets,extent,time,nbin=nbin,width=width,height=height,title=-1,xlabel=xlab,ylabel=ylab,filename=fn,zmax=zmax,timeOffset=timeOffset,cmCustom=cmCustom,isRow=true,keepReg=keepReg,isLast=time==times[end]))
    end
    if savePlot
        save(saveName,g,include_preamble=false)
        return
    else
        return g
    end
end

function plotReachablePD(sets1,sets2,extent;nbin=100,zmax=50,timeOffset=0,cmCustom=true,savePlot=false,saveName="Pics/VCAS_ReachablePD.tex",keepReg=nothing)
    #= Plot a row of reachable plots, which are the reachable regions at different times
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        times (array of int): Times (taus) at which we want to plot (Should be keys in sets)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom viridis colormap
        nbin (int): Number of bins for plotting resolutiones
        saveName (string): Name for the .tex file
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        g (PGFPlots.GroupPlot): Row of reachable region plots
    =#
    PGFPlots.resetPGFPlotsPreamble()
    extraStyle = ",colorbar/draw/.append code={\\begin{axis}[every colorbar,ytick={0},ylabel=,height=.05*\\pgfkeysvalueof{/pgfplots/parent axis width},ymin=0,ymax=0,at={(1.03,0.0)},anchor=south west,colorbar=false]\\end{axis}}"
    g = GroupPlot(2,1, groupStyle = "horizontal sep=0.35cm", style="height=6.0cm, width=6.0cm"*extraStyle)
    width=nothing; height=nothing;
    xlabel=L"$\dot{h}_\text{own}\;(\text{kft/min})$";ylabel=L"$h\;(\text{kft})$"
    
    fn = nothing
    if savePlot
        fn = @sprintf("Pics/VCAS_ReachPD_%d.png",0)
    end
    push!(g,plotReachable(sets1,extent,0,nbin=nbin,width=width,height=height,title=-1,xlabel=xlabel,ylabel=ylabel,filename=fn,zmax=zmax,timeOffset=timeOffset,cmCustom=cmCustom,isRow=true,keepReg=keepReg,isLast=false))
    ylabel=""
    
    fn = nothing
    if savePlot
        fn = @sprintf("Pics/VCAS_ReachPD_%d.png",1)
    end
    push!(g,plotReachable(sets2,extent,0,nbin=nbin,width=width,height=height,title=-1,xlabel=xlabel,ylabel=ylabel,filename=fn,zmax=zmax,timeOffset=timeOffset,cmCustom=cmCustom,isRow=true,keepReg=keepReg,isLast=true))
    ylabel=""
    
    if savePlot
        save(saveName,g,include_preamble=false)
        return
    else
        return g
    end
end

function animateReach(sets,extent,save_name;timeMax=40,timeMin=-40,timeOffset=0,cmCustom=true,zmax=50,nbin=150,keepReg=nothing,time0Num=10)
    #= Create an animation of reachable set plots. Begins with timeMax and counts down to timeMin
    Inputs:
        sets (Dictionary): Reachable sets at different times (taus)
        extent (list of float): Boundaries of plot [xmin, xmax, ymin, ymax]
        timeMin (int): Minimum time (tau) value to plot (Should be a key in  sets)
        timeMax (int): Maximum time (tau) value to plot (Should be a key in  sets)
        save_name (string): Name of gif file (without the .gif)
    Optional Inputs:
        zmax (float): Maximum z-value for colorbar
        timeOffset (int): Offset for time index (tau)
        cmCustom (bool): True if using custom Jet colormap
        nbin (int): Number of bins for plotting resolution
        keepReg (BitArray): Region to always plot as reachable
    Outputs:
        frames (Reel.Frames): Animation of reachable regions changing over time
    =#
    
    currentFrac = 0
    deltaFrac = 0.1
    frames = Frames(MIME("image/svg+xml"), fps=8)
    timeMin = max(timeMin,minimum(keys(sets)))
    
    times = vcat(timeMax:-1:timeMin,zeros(Int64,10).+timeMin)
    if (timeMax>0) .& (timeMin<0)
        times = vcat(timeMax:-1:1,zeros(Int64,time0Num),-1:-1:timeMin,zeros(Int64,10).+timeMin)
    end
    for frame in times
        if (timeMax-frame)/(timeMax-timeMin)>=currentFrac
            println("Completed "*string(round(currentFrac*100))*"%")
            flush(stdout)
            currentFrac+=deltaFrac
        end
        push!(frames,plotReachable(sets,extent,frame,timeOffset=timeOffset,cmCustom=cmCustom,zmax=zmax,nbin=nbin,isLast=frame==timeMin,keepReg=keepReg))
    end
    write(save_name*".gif",frames)
    frames
end
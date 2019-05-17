module Plots

# using Pkg
using ..DISSEQT
using Gadfly
using DataFrames
using Colors

# haskey(Pkg.installed(),"PlotlyJS") && (import PlotlyJS; const py=PlotlyJS) # import to avoid clash with Gadfly


export 
    groupcolors,
    linearregressionplot,
    pairwisescatterplot,
    talusplot,
    fitnesslandscapeplot,
    saveplot


include("colors.jl")
include("stopcodons.jl")
include("pairwise.jl")
include("talus.jl")
include("fitnesslandscape.jl")
include("save.jl")
# haskey(Pkg.installed(),"PlotlyJS") && include("saveplotly.jl")




end
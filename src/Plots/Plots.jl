module Plots

# using Pkg
using ..DISSEQT
using LinearAlgebra
using Gadfly
using DataFrames
using Colors
using Cairo
using Fontconfig

import PlotlyJS; const py=PlotlyJS # import to avoid clash with Gadfly
#using ORCA # Causes lot's of fontconfig warnings on some systems. If user wants to export PlotlyJS figures, then "import ORCA" can be run in the user script insted.


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
include("saveplotly.jl")




end
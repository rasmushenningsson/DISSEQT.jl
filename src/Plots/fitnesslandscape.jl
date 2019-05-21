_plotlycolor(c::AbstractString) = c
_plotlycolor(c) = string("rgb(",round(Int,255*red(c)),',',round(Int,255*green(c)),',',round(Int,255*blue(c)),')')


function logsumexp(x)
    M = maximum(x)
    M + log(sum(exp.(x.-M)))
end


function _alphapoint(model, point::AbstractMatrix, σ::Float64)
    D2 = squareddistances(model.X,point)
    α = logsumexp(-D2/(2σ^2))/log(10)
end
_alpha(model, points::AbstractMatrix, σ::Float64=model.σ) =
    Float64[ _alphapoint(model, points[i:i,:], σ) for i=1:size(points,1) ] # save memory by not computing full D2 matrix

function fitnesslandscapeplot(model, x::AbstractVector{S}, y::AbstractVector{T},
                              points::AbstractMatrix=zeros(0,2), pointDesc=ones(Int,size(points,1)),
                              pointColors=[colorant"black"];
                              width=1024, height=768,
                              zAspect=1,
                              cameraCenter=(0,0,0), cameraEye=(1.25,1.25,1.25),
                              zMin=nothing, zMax=nothing,
                              markerSize=6,
                              σTransparency=0.05) where {S<:Real,T<:Real}
    haskey(Pkg.installed(),"PlotlyJS") || error("Module PlotlyJS not installed, but required for fitnesslandscapeplot().")


    @assert size(points,2)==2
    @assert size(points,1)==length(pointDesc)

    σTransparency *= min(maximum(x)-minimum(x), maximum(y)-minimum(y))

    N,M = length(x), length(y)
    X,Y = repeat(x,1,M), repeat(y',N,1)
    gridPoints = hcat(X[:],Y[:]) # points to evaluate the landscape in
    Z = reshape(predictfitness(model, gridPoints), N, M) # surface
    α = reshape(_alpha(model,gridPoints,σTransparency), N, M) # alpha mask

    if zMin==nothing || zMax==nothing
        zMin,zMax = extrema(Z)
    end


    dz = 0.1*(maximum(Z)-minimum(Z))
    dzP = 0.1dz
    
    # adjust surface depending on α and add white surface hiding the 0-alpha parts
    α1 = clamp.((α+.5)/0.5,0,1)
    α2 = clamp.((α+1.5)/1,0,1)
    Z1 = Z.*α2 + minimum(Z)*(1-α2) - dz*(1-α1)
    Z2 = Z.*α2 + minimum(Z)*(1-α2) - dz*α1

    Z2[[1,end],:],Z2[:,[1,end]] = Z1[[1,end],:],Z1[:,[1,end]] # avoid gaps at the edges

    whiteColorScale = Vector[[0,"rgb(255,255,255)"],[1,"rgb(255,255,255)"]]
    noLighting = py.attr(ambient=1, diffuse=0, specular=0, fresnel=0)

    surf = py.surface(x=X, y=Y, z=Z1, surfacecolor=Z, colorscale="Viridis", cmin=zMin, cmax=zMax, lighting=noLighting)

    surf2 = py.surface(x=X, y=Y, z=Z2, colorscale=whiteColorScale, lighting=noLighting, showscale=false)

    traces = [surf, surf2]

    if !isempty(points)
        pZ = predictfitness(model, points) + dzP

        uniquePointDesc = unique(pointDesc)
        pointIDs = indexin(pointDesc,uniquePointDesc)
        pointColors = String[_plotlycolor(pointColors[i]) for i=uniquePointDesc]
        length(pointColors)==1 && (pointColors = repeat(pointColors,2)) # to avoid degenerate colorscale and div by zero
        nbrColors = length(pointColors)
        colorScale = [ [(i-1)/(nbrColors-1), pointColors[i]] for i=1:nbrColors ]
        col = (pointIDs-1)/(nbrColors-1)

        samples = py.scatter3d(x=points[:,1], y=points[:,2], z=pZ,
                              mode="markers", 
                              marker=py.attr(color=col, colorscale=colorScale, size=markerSize))
        traces = vcat(traces,samples)
    end

    #3-tuples cameraCenter, cameraEye
    camera = py.attr(center=py.attr(x=cameraCenter[1], y=cameraCenter[2], z=cameraCenter[3]),
                     eye=py.attr(x=cameraEye[1], y=cameraEye[2], z=cameraEye[3]))
    scene = py.attr(aspectratio=py.attr(x=1,y=1,z=zAspect),
                    camera=camera, zaxis=py.attr(range=[zMin,zMax]))
    layout = py.Layout(width=width, height=height, scene=scene)

    py.plot(traces, layout)
end




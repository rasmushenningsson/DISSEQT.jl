
function timeseriesindices(t, identifiers...)
    N = length(t)
    @assert all(x->length(x)==N,identifiers)
    combinedIds = collect(zip(identifiers...))
    uniqueIds = unique(combinedIds)
    repInd = indexin(combinedIds, uniqueIds) # assign the same index to all samples in the same replicate

    ind = Vector{Int}[find(repInd.==i) for i=1:length(uniqueIds)] # each entry is a vector with all indices for that replicate
    ind = map(x->x[sortperm(t[x])], ind) # sort each time series by increasing time
    ind, uniqueIds
end

_indexdim(N,dim,ind) = [ i==dim ? ind : Colon() for i=1:N ] # utility to put : as index for every dimension except the chosen one

function timeseries(X, timeSeriesIndices, dim=1)
    @assert dim<=ndims(X)
    map( ind->X[_indexdim(ndims(X),dim,ind)...], timeSeriesIndices )
end



function _interpolate(tStart, tEnd, A, B, t)
    t==tStart && return A
    t==tEnd && return B
    α = (t-tStart)/(tEnd-tStart)
    (1-α)*A + α*B
end


function timeseriesat(keyTimes::AbstractVector{T},V::AbstractMatrix{S},t::Number; clamp=false) where {T,S}
    @assert all(diff(keyTimes).>0)
    if clamp
        t = Base.clamp(t,keyTimes[1],keyTimes[end])
    else
        @assert keyTimes[1] <= t <= keyTimes[end]
    end

    ind = searchsortedlast(keyTimes,t)
    keyTimes[ind]   == t && return V[:,ind]   # avoids unnecessary interpolation and handles edge case
    keyTimes[ind+1] == t && return V[:,ind+1] # avoids unnecessary interpolation
    _interpolate(keyTimes[ind],keyTimes[ind+1], V[:,ind], V[:,ind+1], t)
end

timeseriesat(keyTimes::AbstractVector,V::AbstractMatrix,t::AbstractVector{T}; clamp=false) where {T<:Number} =
    hcat( [timeseriesat(keyTimes,V,s;clamp=clamp) for s in t]... )

timeseriesat(keyTimes::AbstractVector{T},V::AbstractVector{S},i,t; kwargs...) where {T<:AbstractVector,S<:AbstractMatrix} =
    timeseriesat(keyTimes[i],V[i],t; kwargs...)



function timeseriesmean(keyTimes::AbstractVector,V::AbstractMatrix)
    w = diff(keyTimes)/2
    w = vcat(w,0)+vcat(0,w)
    w = w/sum(w)
    sum(V.*w', 2)
end
timeseriesmean(keyTimes::AbstractVector{T},V::AbstractVector{S}) where {T<:AbstractVector,S<:AbstractMatrix} = map(timeseriesmean,keyTimes,V)



# outputs d(V1,V2)(t), the distance of the multivariate time series V1 and V2 (DxT1, DxT2 matrices) as a function of time 
# output given as a piecewise 2nd degree polynomial, 1st row is constant factor, 2nd row is linear factor, 3rd row is square factor (the polynomial describes squared distance)
# each piece is paramterized for the interval [0, tEnd-tStart]
function distancecurve(t1::AbstractVector{T},t2::AbstractVector{T},V1::AbstractMatrix{S},V2::AbstractMatrix{S}) where {T,S}
    #@assert issorted(t1) && issorted(t2)
    @assert all(diff(t1).>0) && all(diff(t2).>0) # check that the values are strictyly increasing

    # do not deal with degenerate cases
    @assert maximum(t1)>minimum(t2) && maximum(t2)>minimum(t1)
    @assert length(t1)>1 && length(t2)>1


    t = Vector{T}()
    params = Vector{S}()

    i=j=1 # index in t1 and t2
    while true
        tStart = max(t1[i], t2[j])
        tEnd = min(t1[i+1], t2[j+1])

        # handle interval
        if tEnd>tStart
            # starting vertices
            A1 = _interpolate(t1[i],t1[i+1],V1[:,i],V1[:,i+1],tStart)
            A2 = _interpolate(t2[j],t2[j+1],V2[:,j],V2[:,j+1],tStart)

            # ending vertices
            B1 = _interpolate(t1[i],t1[i+1],V1[:,i],V1[:,i+1],tEnd)
            B2 = _interpolate(t2[j],t2[j+1],V2[:,j],V2[:,j+1],tEnd)
          
            # l: P+sQ   # s∈[0,tEnd-tStart]
            P = A2-A1   # difference in starting pos
            Q = B2-B1-P # difference in direction
            Q = Q/(tEnd-tStart) # rescale s.t. parameter ranges from 0 to tEnd-tStart

            # ||l(s)||² = ||P+sQ||² = t²||Q||² + 2t<P,Q> + ||P||² = at² + bt + c
            a = vecdot(Q,Q)
            b = 2*vecdot(P,Q)
            c = vecdot(P,P)

            push!(t, tStart)
            push!(params, c, b, a)
        end

        # move one or both of i and j forward
        t1[i+1]==tEnd && (i+=1)
        t2[j+1]==tEnd && (j+=1)

        if i>=length(t1) || j>=length(t2)
            push!(t,tEnd)
            break        
        end
    end

    params = reshape(params, 3, length(t)-1)    
    t, params
end
distancecurve(t1::AbstractVector,t2::AbstractVector,V1::AbstractMatrix,V2::AbstractMatrix) = distancecurve(promote(t1,t2)...,promote(V1,V2)...)



# find all values of t for which f(t)==y, where f(t) is the square root of a positive piecewise 2nd-degree polynomial
function distancecurveequals(t::AbstractVector{T}, P::AbstractVecOrMat, y) where {T}
    @assert y>=0
    y = y^2

    tOut = Vector{promote_type(Float64,T)}()

    P[1,1]==y && push!(tOut, t[1]) # edge case

    for i=1:length(t)-1
        # solve in this segment
        c,b,a = (P[:,i]...,)
        c -= y

        if abs(a)<1e-16
            # linear case;  bx+c=0

            abs(b)<1e-16 && continue # constant case, ignore solutions since would consist of the entire interval

            x = -c/b

            # check that the solution is within the interval
            0<x<=t[i+1]-t[i] && push!(tOut, t[i]+x) 
        else
            # solve 2nd degree polynomial
            w = b^2-4a*c
            w<0 && continue # no solution in this segment

            xmid = -b/(2a)
            xpm  = sqrt(w)/(2a)

            x1 = xmid - xpm
            x2 = xmid + xpm

            # println(i," ",x1," ", x2)

            # check that the solutions are within the interval
            0<x1<=t[i+1]-t[i] &&        push!(tOut, t[i]+x1)
            0<x2<=t[i+1]-t[i] && w>0 && push!(tOut, t[i]+x2) # avoid dubble root for w==0
        end
    end
    tOut
end

# first time the distance reaches a given threshold
function bifurcationtime(t::AbstractVector{T}, P::AbstractVecOrMat, y) where {T}
    RT = promote_type(Float64,T)
    y^2<=P[1,1] && return convert(RT,t[1])
    tEq = distancecurveequals(t,P,y)
    isempty(tEq) ? convert(RT,Inf) : tEq[1] # first time
end



function evaluatedistancecurve(t::AbstractVector, P::AbstractVecOrMat, x::AbstractVector{T}) where {T}
    @assert issorted(t)
    @assert issorted(x)
    @assert size(P,2)==length(t)-1

    y = convert(T,NaN)*ones(x)


    xStartInd = searchsortedfirst(x,t[1])
    

    # for each segment
    for i=1:length(t)-1
        tStart = t[i]
        tEnd = t[i+1]
        xEndInd = searchsortedlast(x,tEnd)

        # evaluate
        if xEndInd>=xStartInd
            xl = x[xStartInd:xEndInd] - tStart # 
            y[xStartInd:xEndInd] = sqrt( [ones(xl) xl xl.^2]*P[:,i] )
        end

        xStartInd = xEndInd+1
    end
    y
end


function distancecurvemean(t::AbstractVector, P::AbstractVecOrMat{T}) where {T}
    # integrate distance

    integral = zero(T)

    # for each segment
    for i=1:length(t)-1
        # find primitive function

        # the 2nd degree polynomial for the segment is by construction nonnegative for all t (not only in the current segment)
        # use formula ∫√(x²+u²)dx = 1/2 x√(x²+u²) + 1/2u²ln(x + √(x²+u²))

        c,b,a = (P[:,i]...,)

        if abs(a)<1e-16
            # linear case
            error("Linear case not yet implemented, but should be very rare.")
        end


        b,c = b/a, c/a # x²+bx+c

        # complete the squares (x+α)²+β
        α = b/2
        β = c-α^2

        β<0 && error("Unexpected negative distance encountered (numerical problems?)")


        # translate function: now √(x²+β)
        # integrate from x0 to x1
        x0 = α
        x1 = α + t[i+1]-t[i]

        sq1 = sqrt(x1^2+β)
        Fx1 = x1*sq1 + β*log(x1+sq1)

        sq0 = sqrt(x0^2+β)
        Fx0 = x0*sq0 + β*log(x0+sq0)

        integral += sqrt(a) * 1/2*(Fx1-Fx0)
    end

    integral/(t[end]-t[1]) # mean
end

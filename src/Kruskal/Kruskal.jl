module Kruskal

using LinearAlgebra
using Statistics
using Optim

export kruskalmds

# naming similar to [Kruskal, Nonmetric Multidimensional Scaling: A Numerical Method]
# X is a vector with values [x₁,y₁,...,x₂,y₂,...,xₙ,yₙ,...], i.e. a vectorized pxN matrix

# gradient-free S, for debugging purposes
# function stressreference{T}(X::AbstractVecOrMat{T},p,D)
#     N = div(length(X),p)

#     num = den = zero(T)

#     # for i>j
#     ind = 1 # index in vectorization of upper triangular part 
#     for j=2:N
#         for i=1:j-1
#             i2 = p*(i-1) # 0-based index of xᵢ in X
#             j2 = p*(j-1) # 0-based index of xⱼ in X

#             t = zero(T)
#             for dim=1:p
#                 t += (X[i2+dim]-X[j2+dim])^2
#             end
#             num += (sqrt(t)-D[ind])^2
#             den += t
            
#             ind += 1
#         end
#     end
#     sqrt(num/den)
# end




mutable struct SData{T}
    lastX::Vector{T} # Vectorization of pxN data point matrix. Needed to figure out whether we moved or not since last call to _update
    p::Int # Number of dimensions in low-dimensional representation.
    N::Int # Number of data points
    targetD::Vector{T} # Vectorization of upper triangular part of target distance matrix.
    num::T
    den::T
    S::T
    numGrad::Vector{T} # gradient of numerator
    denGrad::Vector{T} # gradient of denominator
end
SData(p,N,targetD::Vector{T}) where {T} = SData(NaN*zeros(T,p*N),p,N,targetD,zero(T),zero(T),zero(T),Vector{T}(p*N),Vector{T}(p*N))


function update!(sData::SData{T},X) where {T}
    sData.lastX == X && return # nothing to do

    targetD = sData.targetD
    N = sData.N # number of data points
    p = sData.p
    numGrad = sData.numGrad
    denGrad = sData.denGrad
    num = den = zero(T)
    numGrad[:] = 0
    denGrad[:] = 0
    

    # for i<j, point indices
    ind = 1 # index in vectorization of upper triangular part 
    for j=2:N
        for i=1:j-1
            i2 = p*(i-1) # 0-based index of xᵢ in X
            j2 = p*(j-1) # 0-based index of xⱼ in X

            t = zero(T)
            for dim=1:p
                t += (X[i2+dim]-X[j2+dim])^2
            end
            st = sqrt(t)
            num += (st-targetD[ind])^2
            den += t

            u = 2*(st-targetD[ind])/st # common computations for all dimensions in numGrad
            for dim=1:p
                v = (X[i2+dim]-X[j2+dim])
                #w = u*v
                w = abs(v)<1e-12 ? 0.0 : u*v # handle edge case by letting 0*Inf := 0
                numGrad[i2+dim] += w
                numGrad[j2+dim] -= w
                denGrad[i2+dim] += 2v
                denGrad[j2+dim] -= 2v
            end

            ind += 1
        end
    end


    sData.num = num
    sData.den = den
    sData.S = sqrt(num/den)
    nothing
end


# gradient-aware S (i.e. sharing data with gradient computation)
function stress(X,sData::SData)
    update!(sData,X)
    sData.S
end


function stressgradient(X,storage,sData::SData)
    update!(sData,X)

    targetD = sData.targetD
    num = sData.num
    den = sData.den
    S = sData.S
    numGrad = sData.numGrad
    denGrad = sData.denGrad

    C = 1.0/(2*S) # constant factor
    for i=1:length(X)
        storage[i] = C*(numGrad[i]*den - num*denGrad[i]) / den^2
    end

    nothing
end




function updatelocation(D::AbstractVector{T}, X0::AbstractMatrix{T}, tol) where {T}
    p = size(X0,1)
    N = size(X0,2)
    X0 = X0[:]

    # # Simple, gradient-free version
    # res = optimize(x->stressreference(x,p,D), X0, BFGS())
    # @assert Optim.converged(res)
    # X = reshape(Optim.minimizer(res),p,N)


    sData = SData(p,N,D) # init storage for reused computations
    # res = optimize(x->stress(x,sData), (storage,x)->stressgradient(x,storage,sData), X0, BFGS())
    res = optimize(x->stress(x,sData), (storage,x)->stressgradient(x,storage,sData), X0, BFGS(), Optim.Options(x_tol=tol, f_tol=tol))
    #res = optimize(x->stress(x,sData), X0, BFGS()) # doesn't use gradient, for comparison
    Optim.converged(res) || error("Kruskal MDS: Location optimization did not coverge.")
    X = reshape(Optim.minimizer(res),p,N)


    # Since the optimum is preserved when X is uniformly scaled, normalize scaling to keep numerical stability
    X = X .- mean(X;dims=2) # center X (doesn't change stress)
    X = X / sqrt(mean(sum(X.^2;dims=1);dims=2)[1]) # scale so that RMS distance to origin is equal to 1

    X, Optim.minimum(res)
end



# convert distance matrix to pdist-style vector
function distancevector(D::AbstractMatrix)
    # vectorize upper triangular part of D
    N = size(D,1)
    dist = zeros(div(N*(N-1),2))
    ind=1
    for j=2:N
        for i=1:j-1
            dist[ind] = D[i,j]
            ind += 1
        end
    end
    dist
end


# convert pdist-style distance vector to matrix
function distancematrix!(D::AbstractMatrix, dist::AbstractVector)
    N = size(D,1)

    for i=1:N
        D[i,i] = 0.0
    end

    ind=1
    for j=2:N
        for i=1:j-1
            D[i,j] = D[j,i] = dist[ind]
            ind += 1
        end
    end
    D
end
distancematrix(dist::AbstractVector{T}) where {T} = (M=length(dist); N=round(Int,1/2+sqrt(1/4+2M)); distancematrix!(Matrix{T}(N,N),dist))



function isotonicregression(Y::AbstractVector{T}, tieGroups::Vector{UnitRange}=[]) where {T}
    X = Tuple{T,Int}[] # tuple with values and integer weight

    # reorder values in each tieGroup to put them in increasing order
    isempty(tieGroups) || (Y=copy(Y)) # avoid modifying the original vector
    tiePerms = map!(rng->sortperm(Y[rng]), Vector{Vector{Int}}(undef,length(tieGroups)), tieGroups)
    for (g,p) in zip(tieGroups,tiePerms)
        Y[g] = Y[g[p]]
    end


    for y in Y
        yw = 1

        while !isempty(X) && y<=X[end][1]
            # merge with previous point
            x,xw = pop!(X)
            y,yw = (x*xw+y*yw)/(xw+yw), xw+yw
        end
        push!(X, (y,yw)) # store point
    end

    # unwrap to full vector by duplicating elements by their weight
    Z = Vector{T}(undef,length(Y))
    zi = 1
    for x in X
        Z[zi:zi+x[2]-1] .= x[1]
        zi += x[2]
    end


    # inverse tieGroup reordering
    for (g,p) in zip(tieGroups,tiePerms)
        Z[g[p]] = Z[g]
    end

    Z
end


function pdist(X)
    # create pdist-style distance vector
    p = size(X,1)
    N = size(X,2)
    dist = similar(X, div(N*(N-1),2))
    ind=1
    for j=2:N
        for i=1:j-1
            dist[ind] = sqrt(sum((X[:,i]-X[:,j]).^2)) # D[i,j]
            ind += 1
        end
    end
    dist
end


function duplicateranges(v)
    @assert issorted(v)
    count = 0
    out = UnitRange[]
    for i=1:length(v)
        if i<length(v) && v[i]==v[i+1]
            count+=1
        else
            count>0 && push!(out, i-count:i)
            count = 0
        end
    end
    out
end


function kruskalmds(D::AbstractMatrix{T}, p::Integer, X0=convert(Matrix{T},randn(size(D,1),p)); maxIter=1000, tol=1e-6) where {T<:AbstractFloat}
    @assert p>0
    @assert D≈D' "Distance matrix must be symmetric"
    N = size(D,1)
    @assert size(X0)==(N,p)
    X0 = X0' # work with transpose internally.

    dist = distancevector(D)
    perm = sortperm(dist)
    iperm = invperm(perm)
    tieGroups = duplicateranges(dist[perm])


    X = X0 .- mean(X0;dims=2) # center X (doesn't change stress)
    X = X / sqrt(mean(sum(X.^2;dims=1);dims=2)[1]) # scale so that RMS distance to origin is equal to 1 (doesn't change stress, if the monotone function is updated accordingly)


    S = Inf
    converged = false


    # println("i: 0, S=", Sreference(X0,p,dist))

    for i=1:maxIter # nbr iterations
        prevS = S
        modelDist = pdist(X)[perm] # compute distance and permute
        permDist = isotonicregression(modelDist, tieGroups) # nonnegativity constraint will be enforced automatically since all values in modelDist are nonnegative
        dist = permDist[iperm]

        XOld = X
        X,S = updatelocation(dist,X,tol)
        stepSize = vecnorm(X-XOld)

        # println("i: ", i, ", S=", S, ", stepSize=", stepSize, ", prevS-S=", prevS-S, ", tol*S=", tol*S)
        stepSize<tol && (converged=true; break)
        prevS-S<tol*S && (converged=true; break)
        prevS = S
    end

    # post-processing
    # the stress is invariant under rotation, so we can normalize rotation by principal axes
    F = svd(X)
    X = F.V*Diagonal(F.S)

    X, S, converged, distancematrix(dist)
end
kruskalmds(D::AbstractMatrix, p::Integer, X0=randn(size(D,1),p); kwargs...) = kruskalmds(convert(Matrix{Float64},D),p,X0;kwargs...)

end

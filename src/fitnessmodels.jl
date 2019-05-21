

# returns MxN matrix of squared distance values, where M is the number of test points and N the number of model points
function squareddistances(X::Matrix{Float64}, Y::Matrix{Float64})
    @assert size(Y,2)==size(X,2) # sequence space dimensions should be the same
    sum(Y.^2;dims=2)' .+ sum(X.^2;dims=2) .- 2*X*Y'
end

# returns MxN matrix of edit distance values
function editdistance(X,Y)
    @assert size(Y,2)==size(X,2)
    M,N = size(X,1),size(Y,1)
    D = zeros(Int,M,N)
    for j=1:N
        for i=1:M
            D[i,j] = sum(X[i,:].!=Y[j,:])
        end
    end
    D
end

function groupequal(X,Y)
    @assert size(Y,2)==size(X,2)
    M,N = size(X,1),size(Y,1)
    D = BitArray(M,N)
    for j=1:N
        for i=1:M
            D[i,j] = X[i,:] == Y[j,:]
        end
    end
    D
end


mutable struct LandscapeModel
    X::Matrix{Float64}
    fitness::Vector{Float64}
    σ::Float64

    function LandscapeModel(X::Matrix{Float64}, fitness::Vector{Float64}, σ::Float64)
        @assert size(X,1)==length(fitness)
        new(X,fitness,σ)
    end
end
LandscapeModel(X::AbstractMatrix,fitness::AbstractVector,σ::Real) = LandscapeModel(convert(Matrix{Float64},X), convert(Vector{Float64},fitness), convert(Float64,σ))

LandscapeModel(X::Matrix{Float64},fitness::Vector{Float64},σValues::AbstractVector{Float64}; kwargs...) = 
    LandscapeModel(X,fitness,landscapekernelwidth(X,fitness,σValues;kwargs...)[1])
LandscapeModel(X,fitness,σValues::AbstractVector; kwargs...) =
    LandscapeModel(convert(Matrix{Float64},X), convert(Vector{Float64},fitness), convert(Vector{Float64},σValues); kwargs...)





# D2 - MxN matrix of squared distance values, where M is the number of test points and N the number of training points
function _predictfitness(D2::Matrix{Float64}, fitness::Vector{Float64}, σ::Float64)
    # Wᵢⱼ = e^( -Dᵢⱼ²/(2σ²) )
    # But since we only care about each row of W up to a constant value, we can multiply each row with e^(Aᵢ/2σ²) where Aᵢ=minⱼ Dᵢⱼ²
    # this guarantees that we get precision in the computations, by making Wᵢⱼ=1 for the closest training point.
    W = exp.( (minimum(D2;dims=2).-D2)/(2σ^2) )
    W = W./sum(W;dims=2) # normalize sum of weights to 1 for each row
    W*fitness # the test point fitness is a weighted average of the fitness values for the training points
end

predictfitness(model::LandscapeModel, X::Matrix{Float64}) = _predictfitness(squareddistances(X,model.X), model.fitness, model.σ)



_fitnessresidualssquared(D2::Matrix{Float64}, fitnessTrain::Vector{Float64}, fitnessTest::Vector{Float64}, σ::Float64) =
    (fitnessTest - _predictfitness(D2, fitnessTrain, σ)).^2


_fitnesserrortransform(x) = x       # mean squared error
#_fitnesserrortransform(x) = sqrt(x) # RMS error

_fitnesserror(residualsSquared::Array{Float64}, args...; kwargs...) = _fitnesserrortransform(mean(residualsSquared, args...; kwargs...))

_fitnesserror(D2::Matrix{Float64}, fitnessTrain::Vector{Float64}, fitnessTest::Vector{Float64}, σ::Float64) =
    _fitnesserror(_fitnessresidualssquared(D2,fitnessTrain,fitnessTest,σ))


# find optimal kernel width by repeated random subsampling
function landscapekernelwidth(X::Matrix{Float64}, fitness::Vector{Float64}, σValues::AbstractVector{Float64}; nbrIter=100, trainingProp=0.9)
    N = size(X,1)
    X2 = sum(X.^2; dims=2)
    D2 = X2 .+ X2' .- 2*X*X'

    nbrTest = round(Int,N*(1-trainingProp))
    @assert 0 < nbrTest < N

    allErrors = zeros(length(σValues), nbrIter)

    testMask = falses(N)
    for i=1:nbrIter
        # select a random training subset each time
        testMask[:] .= false
        testMask[sample(1:N,nbrTest,replace=false)] .= true

        D2TestTrain = D2[testMask,.~testMask]
        fitnessTrain = fitness[.~testMask]
        fitnessTest = fitness[testMask]

        allErrors[:,i] = map(x->_fitnesserror(D2TestTrain,fitnessTrain,fitnessTest,x), σValues)
    end

    meanErrors = mean(allErrors,dims=2)[:]
    ind = argmin(meanErrors)
    σValues[ind], meanErrors[ind], meanErrors
end





# find optimal kernel width by repeated random subsampling
# gives one σ per sample, where σ is estimated from all samples except the one left out
function leaveoneoutlandscapekernelwidth(X::Matrix{Float64}, fitness::Vector{Float64}, σValues::AbstractVector{Float64}; nbrIter=100, trainingProp=0.9)
    N = size(X,1)
    @assert N==length(fitness)
    X2 = sum(X.^2;dims=2)
    D2 = X2 .+ X2' .- 2*X*X'

    nbrTest = round(Int,N*(1-trainingProp))
    @assert 1 < nbrTest < N

    #allErrors = zeros(length(σValues), nbrIter)
    allErrors = zeros(N, length(σValues), nbrIter)

    testMask = falses(N)
    for i=1:nbrIter
        # select a random training subset each time
        testMask[:] = false
        testMask[sample(1:N,nbrTest,replace=false)] = true

        D2TestTrain = D2[testMask,.~testMask]
        fitnessTrain = fitness[.~testMask]
        fitnessTest = fitness[testMask]


        # TODO: merge common computations of case A and B


        # Case A, samples part of the Training set (different set of weights depending on which sample we leave out)
        for (j,σ) in enumerate(σValues)
            W = exp.( (minimum(D2TestTrain;dims=2).-D2TestTrain)/(2σ^2) )

            A = W*fitnessTrain .- W.*fitnessTrain' # nbrTest x nbrTrain - unnormalized fitness for test sample i, with training sample j removed
            B = sum(W;dims=2) .- W                 # nbrTest x nbrTrain - sum of weights for test sample i, with training sample j removed

            # for elements with Bᵤᵥ<<1, we are not guaranteed good precision in the computations (because the closest sample has a much larger weight than the others, and it has been removed)
            # update Aᵤᵥ and Bᵤᵥ for those elements
            # for (u,v) in zip(ind2sub(size(B), findall(B.<1e-3))...) # 1e-3 is still a conservative bound
            for (u,v) in zip(Tuple.(CartesianIndices(B)[findall(B.<1e-3)])...) # 1e-3 is still a conservative bound
                d2 = D2TestTrain[u,:]
                d2[v] = Inf # remove training sample by setting it to infinitely far away
                w = exp.( (minimum(d2).-d2)/(2σ^2) ) # note that minimum(d2) is different from above, which gives us the precision needed
                A[u,v] = dot(w,fitnessTrain)
                B[u,v] = sum(w)
            end

            #minimum(B)<0.99 && println(minimum(B)) # printout when the sample with the largest weight has been removed, but precision was good enough
            #@assert all(B.>4eps()) # just for testing the code

            F = A./B # nbrTest x nbrTrain - Fᵤᵥ = predicted value of test sample i, with training sample j removed
            residuals = (F .- fitnessTest).^2 # nbrTest x nbrTrain - squared residual for test sample i, with training sample j removed

            #allErrors[.~testMask, j, i] = dropdims(mean(residuals;dims=1); dims=1) # nbrTrain - mean squared error with training sample j removed
            allErrors[.~testMask, j, i] = dropdims(_fitnesserror(residuals;dims=1); dims=1) # nbrTrain - mean squared error with training sample j removed
        end
        
        # Case B, samples part of the Test set (different set of residuals depending on which sample we leave out)
        for (j,σ) in enumerate(σValues)
            residuals = _fitnessresidualssquared(D2TestTrain, fitnessTrain, fitnessTest, σ)
            totError = sum(residuals)

            # take the mean with one sample removed
            # allErrors[testMask, j, i] = (totError.-residuals)/(nbrTest-1)
            allErrors[testMask, j, i] = _fitnesserrortransform( (totError.-residuals)/(nbrTest-1) )
        end
    end

    meanErrors = dropdims(mean(allErrors; dims=3); dims=3)
    ind = mapslices(argmin, meanErrors; dims=2)[:]
    σValues[ind], getindex.((meanErrors,), 1:length(ind), ind), meanErrors
end









mutable struct NearestNeighborModel{T}
    X::Matrix{T}
    fitness::Vector{Float64}
    function NearestNeighborModel{T}(X::Matrix{T}, fitness::Vector{Float64}) where T
        @assert size(X,1)==length(fitness)
        new(X,fitness)
    end
end
# euclidean distance
NearestNeighborModel(X::AbstractMatrix{T},fitness) where {T<:Number} = NearestNeighborModel{Float64}(convert(Matrix{Float64},X), convert(Vector{Float64},fitness))
# edit distance
NearestNeighborModel(X::AbstractMatrix{Symbol},fitness) = NearestNeighborModel{Symbol}(convert(Matrix{Symbol},X), convert(Vector{Float64},fitness))

# euclidean distance
_predictfitnessnn(D2::Matrix{Float64}, fitness::Vector{Float64}) =
    Float64[ fitness[argmin(D2[i,:])] for i=1:size(D2,1) ] # find nearest neighbour (expect no ties for Float64)
predictfitness(model::NearestNeighborModel{Float64}, X::Matrix{Float64}) = _predictfitnessnn(squareddistances(X,model.X), model.fitness)

# edit distance
function _predictfitnesseditdist(D::Matrix{Int}, fitness::Vector{Float64})
    M,N = size(D)
    minValues = minimum(D;dims=2)
    hits = D.==minValues # mask of all closest samples
    Float64[ mean(fitness[hits[i,:][:]]) for i=1:M ]
end
predictfitness(model::NearestNeighborModel{Symbol}, X::Matrix{Symbol}) = _predictfitnesseditdist(editdistance(X,model.X), model.fitness)




# every row of X should consist of categorical values will be compared to other rows
# X could be for instance be a vector of group IDs, a Matrix with categorical values or a DataFrame.
mutable struct GroupModel{T}
    X::T
    fitness::Vector{Float64}
    function GroupModel{T}(X::T, fitness::Vector{Float64}) where T
        @assert size(X,1)==length(fitness)
        new(X,fitness)
    end
end
GroupModel(X::T,fitness) where {T} = GroupModel{T}(X, convert(Vector{Float64},fitness))

function _predictfitnessgroup(D::AbstractArray{Bool}, fitness::Vector{Float64})
    M,N = size(D)
    f = fill(NaN,M)
    for i=1:M
        any(D[i,:]) || continue # return NaN if there is no matching group
        f[i] = mean( fitness[D[i,:][:]] )
    end
    f
end
predictfitness(model::GroupModel{T}, X::T) where {T} = _predictfitnessgroup(groupequal(X,model.X), model.fitness)




function leaveoneoutpredict(modelType, modelData, fitness, modelArgs...; perSampleModelArgs=())
    typeof(perSampleModelArgs) <: Tuple || (perSampleModelArgs=(perSampleModelArgs,)) # wrap in tuple if a single argument was given

    N = size(modelData,1)
    nd = ndims(modelData)
    colons = repeat([Colon()],nd-1) # for slicing in first dimension only

    predicted = zeros(N)
    for i=1:N
        mask = trues(N)
        mask[i] = false # all except the current one
        model = modelType(modelData[mask,colons...], fitness[mask], modelArgs..., map(x->x[i], perSampleModelArgs)...)
        predicted[i] = predictfitness(model, modelData[i:i,colons...])[1]
    end
    predicted
end


varianceunexplained(predicted,truth) = sum((predicted-truth).^2) / sum((truth-mean(truth)).^2)

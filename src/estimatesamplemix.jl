
_uniqueid(a::AbstractVector) = indexin(a,unique(a)) # TODO: replace with indunique() or similar when it is introduced in julia (c.f. sort() and sortperm())

function _partitionpositions(referenceCodons::AbstractArray{Symbol,2})
    # for each position, assign an ID number to each reference s.t. they have the same ID if they coincide at that position
    IDPerPos = Vector{Int}[ _uniqueid(referenceCodons[i,:][:]) for i=1:size(referenceCodons,1) ]
    _uniqueid(IDPerPos) # give each position a number depending on which equivalence class it is in
end



function _partitionmatrices(partitionInd::AbstractVector{Int}, referenceCodons::AbstractArray{Symbol,2}, freqs::AbstractArray{Float64,3}, codons::Vector{Symbol})
    @assert size(referenceCodons,1)==size(freqs,2)

    nbrSamples = size(freqs,3)
    nbrReferences = size(referenceCodons,2)
    nbrPositions = length(partitionInd)


    nbrUniqueCodons = length(unique(referenceCodons[partitionInd[1],:])) # the number of unique codons is the same at all positions within a partition
    nbrUniqueCodons==1 && return zeros(0,nbrReferences), zeros(0,nbrSamples) # skip partition if it only has one unique codon per site (i.e. the partition doesn't help to determine the sample mix)


    f = zeros(nbrUniqueCodons,nbrPositions,nbrSamples) # UniqueCodons x Positions x Samples
    for (i,p) in enumerate(partitionInd) # for each position in the partition
        c = unique(referenceCodons[p,:]) # the codons belonging to the partition at the given site
        cInd = indexin(c, codons)
        f[:,i,:] = freqs[cInd,p,:] # extract slice for UniqueCodons x Samples
    end

    
    f = f ./ sum(f,1) # normalize freqs to sum to one for each position and sample (i.e. conditional probability given that we have to chose one of the reference codons)


    # output trimmed frequencies (Y matrix)
    Y = zeros(nbrUniqueCodons,nbrSamples)
    for s=1:nbrSamples # for each sample
        # create masks for trimming
        tMask = trues(nbrPositions)
        for k=1:nbrUniqueCodons
            fs = f[k,:,s] # freqs associated to each unqiue group of reference genomes and sample
            q = quantile(fs, [0.25,0.75])
            tMask .&= tMask .& (q[1].<=fs) .& (fs.<=q[2])
        end
        # compute trimmed means
        Y[:,s] = mean(f[:,tMask,s], 2)
    end


    # output reference->codon identification matrix (A matrix)
    # A_ij == 1 iff reference j has codon i and 0 otherwise
    A = zeros(nbrUniqueCodons,nbrReferences)
    A[sub2ind(size(A),_uniqueid(referenceCodons[partitionInd[1],:][:]),1:nbrReferences)] = 1 # figure out which references map to the same codon

    # weight by number of positions in this group
    A *= nbrPositions
    Y *= nbrPositions

    A, Y
end




# TODO: support referenceGenomes::Array{Reference} too
function estimatemix(swarms::AnnotatedArray, referenceGenomes::Vector{T}) where {T<:AbstractString}
    @assert all( map(length,referenceGenomes) .== length(referenceGenomes[1]) )
    
    nbrReferences = length(referenceGenomes)
    nbrSamples = size(swarms,5)

    pos = swarms[:position][:]
    referenceCodons = hcat(seq2codons(referenceGenomes)...) # assumes all references have equal length
    referenceCodons = referenceCodons[pos,:]


    codons = reshape(swarms[:codon],64)
    freqs = reshape(convert(Array, swarms), 64, size(swarms,4), nbrSamples) # Codons x Positions x Samples


    partitions = _partitionpositions(referenceCodons)
    nbrPartitions = maximum(partitions)

    A = zeros(0,nbrReferences)
    Y = zeros(0,nbrSamples)

    for partition=1:nbrPartitions
        Ap, Yp = _partitionmatrices(find(partition.==partitions), referenceCodons, freqs, codons)
        # concat
        A = vcat(A,Ap)
        Y = vcat(Y,Yp)
    end

    # solve separately for each sample for increased precision
    mixes = zeros(nbrSamples,nbrReferences)
    X = Variable(nbrReferences)
    for s=1:nbrSamples
        problem = minimize( sumsquares(A*X-Y[:,s]), [X>=0, ones(nbrReferences)'X==1] )
        solve!(problem, SCSSolver(eps=1e-6,verbose=false), verbose=false)
        problem.status==:Optimal || error("Failed to compute mix for sample ", swarms[:id][s])
        mixes[s,:] = (z=max(X.value,0); z/sum(z)) # the solver might return slightly negative values, get rid of those
    end

    # # Fixing Y should improve performance, but the code below doesn't run faster even when we have hundreds of samples...
    # # solve separately for each sample for increased precision
    # mixes = zeros(nbrSamples,nbrReferences)
    # X = Variable(nbrReferences)
    # y = Variable(size(A,1))
    # problem = minimize( sumsquares(A*X-y), [X>=0, ones(nbrReferences)'X==1] )
    # for s=1:nbrSamples
    #     fix!(y,Y[:,s])
    #     solve!(problem, SCSSolver(eps=1e-6,verbose=false), verbose=false)#, warmstart=s>1)
    #     problem.status==:Optimal || error("Failed to compute mix for sample ", swarms[:id][s])
    #     mixes[s,:] = max(X.value,0) # the solver might return slightly negative values, get rid of those
    # end

    mixes
end

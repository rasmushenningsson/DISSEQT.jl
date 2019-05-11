

function _minimum!(D,clusterInd)
    N = size(D,1)
    m = Inf
    mI = mJ = -1
    for j=2:N
        clusterInd[j]==0 && continue
        for i=1:j-1
            clusterInd[i]==0 && continue
            if D[i,j]<=m
                m = D[i,j]
                mI,mJ = i,j
            end
        end
    end
    m, mI, mJ
end



# only upper triangular part of D is used 
function hierarchicalclustering(D::AbstractMatrix{T}; similarity=false) where {T}
    @assert size(D,1)==size(D,2)
    N = size(D,1)

    distType = promote_type(Float64,T)
    D = copy(convert(Array{distType},D)) # since D will be updated when clusters are merged
    similarity && (D*=-1) # just change sign when using similarity instead of dissimilarity! (Nothing breaks if we have distances below 0. Relative values are all that matter.)


    # index, index, nbr elements in cluster, distance/similarity
    # indices 1-N
    clusters = Vector{Tuple{Int,Int,Int,distType}}()
    
    # the cluster index for each row/column in D. Indices from 1:N are leaves.
    # From N+1 and upwards are in clusters Vector (subtract N from index).
    # A value of 0 means that the row/column is no longer used. (One row/column is discarded every time two clusters are merged.)
    clusterInd = collect(1:N)


    for k=1:N-1 # each created cluster has two children, i.e. we need N-1 clusters (not counting leaves)
        m,mI,mJ = _minimum!(D,clusterInd)

        # create new cluster
        cI,cJ = clusterInd[mI], clusterInd[mJ]

        nbrElementsI = cI<=N ? 1 : clusters[cI-N][3]
        nbrElementsJ = cJ<=N ? 1 : clusters[cJ-N][3]

        dist = D[mI,mJ]
        similarity && (dist*=-1) # change back from distance to similarity if needed
        push!(clusters, (cI,cJ,nbrElementsI+nbrElementsJ,dist))

        # update distance matrix and clusterInd
        clusterInd[mJ] = 0
        clusterInd[mI] = k+N


        # construct new distances from D[mI,u], D[mJ,u]
        for u=1:N
            u==mI && continue
            clusterInd[u]==0 && continue # includes u==mJ

            mI2,u2 = minmax(mI,u) # upper triangular index
            D[mI2,u2] = (D[mI2,u2]*nbrElementsI + D[minmax(mJ,u)...]*nbrElementsJ) / (nbrElementsI+nbrElementsJ) # TODO: change here for different linkage criteria
        end

    end

    clusters
end



function samplesinclusters(clusters)
    N = length(clusters)+1 # Since every merge makes two clusters into one, we can figure out the original number of samples this way.
    samples = Vector{Vector{Int}}(length(clusters))
    for (i,c) in enumerate(clusters)
        samples[i] = vcat( c[1]>N ? samples[c[1]-N] : [c[1]], 
                           c[2]>N ? samples[c[2]-N] : [c[2]] )
    end
    samples
end




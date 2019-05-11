

# get "super node" of a node and update all nodes in the path to the supernode
# Definition: A "super node" is a node in the same connected component, with a lower index than the given node (or equal if the index is the lowest one in the connected component)
function supernode!(superNodes, i)
    if i != superNodes[i]
        superNodes[i] = supernode!(superNodes, superNodes[i])
    end
    superNodes[i]
end


# makes matrix connected
# assumes D and DFull are symmetric
function makeconnected!(D::AbstractMatrix{T},DFull::AbstractMatrix{T}) where {T<:Real}
    N = size(D,1)

    # find all connected components
    superNodes = collect(1:N) # a "super node" is another node in the same connected component, s.t. if we follow through a path of super nodes, we will always reach the same "root" vertex
    for j=2:N
        for i=1:j-1
            D[i,j] != typemax(T) && (superNodes[j] = min(superNodes[j],i)) # always choose the lower one as super node by making the higher index refer to the lower
        end
    end

    # make sure all connected components are represented by their super nodes
    for i=2:N
        superNodes[i] = superNodes[superNodes[i]] # makes sure every index up to i has the lowest index in the connected component as their supernode
    end
    nbrComponents = length(unique(superNodes))

    if nbrComponents>1
        # collect and sort edges
        M = div(N*(N-1),2)
        E    = Vector{T}(M)
        from = Vector{Int}(M)
        to   = Vector{Int}(M)
        ind = 1
        for j=2:N
            for i=1:j-1
                E[ind]    = DFull[i,j]
                from[ind] = i
                to[ind]   = j
                ind += 1
            end
        end
        perm = sortperm(E)
        E    = E[perm]
        from = from[perm]
        to   = to[perm]

        # add edges from shortest to longest, but only if they connect two different components
        for k=1:M
            nbrComponents==1 && break # done

            i = supernode!(superNodes,from[k])
            j = supernode!(superNodes,to[k])
            i==j && continue # the nodes belong to the same component
            i,j = minmax(i,j) 

            # use this edge to connect the components
            superNodes[j] = i
            D[from[k],to[k]] = D[to[k],from[k]] = E[k]
            nbrComponents -= 1
        end
    end
    D
end




# Creates distance matrix based on Isomap distances
function isomapmatrix(D::AbstractMatrix{T}, nearestNeighbours::Integer, distanceThreshold::Real) where {T<:Real}
    @assert ndims(D)==2
    @assert size(D,1)==size(D,2)
    @assert nearestNeighbours>0
    N = size(D,1)
    nearestNeighbours = min(N,nearestNeighbours+1) # add 1 to include self as neighbour

    I = similar(D) # Isomap matrix
    I[:] = typemax(T)

    for v=1:N # for each vertex, look at outgoing edges (i.e. row), and keep all edges in the neighbourhood (i.e. below distanceThreshold, but at least nearestNeighbour edges)
        # distsSorted = sort(D[v,:])
        # k = maximum(nearestNeighbours+1, searchsortedlast(distsSorted,distanceThreshold)) # +1 to include 0-distance to the vertex itself
        # d = distsSorted[k]

        # mask = D[v,:].<=d
        # I[v,mask] = D[v,mask]

        mask = D[v,:] .<= distanceThreshold
        
        if countnz(mask)<nearestNeighbours # not enough neighbours within radius? Include at least nearestNeighbour edges.
            mask = D[v,:] .<= sort(D[v,:][:])[nearestNeighbours] # this includes ties
        end

        I[v,mask] = D[v,mask] # copy edges
    end

    # make symmetric
    I = min.(I,I')


    # make sure graph is connected by adding shortest distance between components
    makeconnected!(I,D)

    # make the graph complete by finding shortest distances between all nodes
    floydwarshall!(I)
end






function isomap(D::AbstractMatrix{T}, p::Integer, nearestNeighbours::Integer, distanceThreshold::Real) where {T<:Real}
    I = isomapmatrix(D, nearestNeighbours, distanceThreshold)
    mds(I,p)
end

function kruskalisomap(D::AbstractMatrix{T}, p::Integer, nearestNeighbours::Integer, distanceThreshold::Real; kwargs...) where {T<:Real}
    I = isomapmatrix(D, nearestNeighbours, distanceThreshold)
    X0 = mds(I,p)
    kruskalmds(I, p, X0; kwargs...)
end
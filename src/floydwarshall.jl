

# saturated addition for positive numbers
function addsaturated{T}(a::T,b::T)
    r = a+b
    r<a || r<b ? typemax(T) : r # check for overflow
end
addsaturated{T<:AbstractFloat}(a::T,b::T) = a+b # Inf+a = Inf, so no overflow check is needed

# Finds the shortest paths between all pairs of nodes in an directed graph, specified by a matrix of pairwise distances.
# If there are no edges between two nodes, the distance should be typemax(T) (which is Inf for Float64).
function floydwarshall!{T<:Real}(D::AbstractMatrix{T})
    @assert ndims(D)==2
    @assert size(D,1)==size(D,2)
    @assert all(D.>=0)
    N = size(D,1)


    for k=1:N
        for j=1:N
            for i=1:N
                #D[j,i] = min(D[j,i], D[j,k]+D[k,i])
                D[j,i] = min(D[j,i], addsaturated(D[j,k],D[k,i]))
            end
        end
    end

    # This is not needed since we assume the matrix has only positive entries
    # for i=1:N
    #     @assert D[i,i]>=0 "The graph has negative cycles."
    # end
    D
end
floydwarshall{T<:Real}(D::AbstractMatrix{T}) = floydwarshall!(copy(D))



mutable struct MedianAccumulatorSimple{T}
    values::Vector{T}
    ind::Int
end
MedianAccumulatorSimple(::Type{T},N) where {T} = MedianAccumulatorSimple(Vector{T}(2N+1),1)


function init!(acc::MedianAccumulatorSimple{T}, A::AbstractArray{T}, N::Integer, dim::Integer, Ipost::CartesianIndex, Ipre::CartesianIndex) where {T}
    D = size(A,dim)

    ind = 1
    for i=1-N:1 # N elements before and the starting element
        acc.values[ind] = A[Ipre,1,Ipost] # duplicate initial edge element
        ind+=1
    end
    for i=2:N+1 # N elements after the starting element
        acc.values[ind] = A[Ipre,min(i,D),Ipost] # duplicate edge elements (needed for large N)
        ind += 1
    end
    acc.ind = 1
    nothing
end
function update!(acc::MedianAccumulatorSimple{T}, new::T) where {T}
    acc.values[acc.ind] = new
    acc.ind = acc.ind==length(acc.values) ? 1 : acc.ind+1
    nothing
end
_median(acc::MedianAccumulatorSimple) = median(acc.values)


# 1-dimensional median filter along the selected dimension
# elements outside the range are considered to be the same as the edge elements
function medianfilter(A::AbstractArray{T}, N::Integer, dim::Integer=1) where {T}
    @assert N>=0
    (N == 0 || size(A,dim)==0) && return copy(A) # no-op
    B = similar(A)

    # Cartesian indexing follows exponential filter example from http://www.juliabloggers.com/multidimensional-algorithms-and-iteration/ 
    Rpre = CartesianRange(size(A)[1:dim-1])
    Rpost = CartesianRange(size(A)[dim+1:end])

    D = size(A,dim)

    acc = MedianAccumulatorSimple(T,N) # possibly have one accumulator for everything in Rpre to improve cache efficiency?

    for Ipost in Rpost
        for Ipre in Rpre
            # initialize
            init!(acc, A, N, dim, Ipost, Ipre)

            # moving filter
            i=1
            while true
                i1 = i-N
                i2 = i+N+1
                B[Ipre,i,Ipost] = _median(acc)

                i == size(A,dim) && break
                # update!(acc, A[Ipre,max(i1,1),Ipost], i1, A[Ipre,min(i2,D),Ipost], i2)
                update!(acc, A[Ipre,min(i2,D),Ipost])
                i += 1
            end
        end
    end

    B
end


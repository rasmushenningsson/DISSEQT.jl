

# assumes symmetric matrix
function mds{T<:AbstractFloat}(D::AbstractMatrix{T}, p::Integer)
    @assert D==D'
    N = size(D,1)
    
    D = D.^2 # MDS works with squared distances

    m = mean(D,1)
    D = Symmetric(-1/2*(D.-m.-m'+mean(m))) # remove mean

    Σ,X = eig(D, N-p+1:N)

    #X = diagm(sqrt(Σ))*X'
    #X[p:-1:1,:] # reorder since eig orders from smallest to largest eigenvalue and we want largest to smallest
    X = X*diagm(sqrt.(Σ))
    X[:,p:-1:1] # reorder since eig orders from smallest to largest eigenvalue and we want largest to smallest
end

mds(D::AbstractMatrix, p::Integer) = mds(convert(Matrix{Float64},D),p)


function talusplot(σ2::Vector, k=1:length(σ2)-1, args...; kwargs...)
    talus = -diff(log2.(σ2))/2 # log₂(σₖ)-log₂(σₖ₊₁)

    df = DataFrame(k=k, talus=talus[k])

    plot(layer(df, x=:k, y=:talus, Geom.line()),
         layer(df, x=:k, y=:talus, Geom.point()),
         Guide.xlabel("k (Principal Component)"), Guide.ylabel("log₂(σₖ)-log₂(σₖ₊₁)"),
         args...; kwargs...)
end


function talusplot(X::Matrix, k=[], args...; kwargs...)
    K = size(X,2)<=size(X,1) ? Symmetric(X'X) : Symmetric(X*X') # Work with the smaller of the kernel matrices

    # compute σᵢ² and talus plot
    σ2 = eigvals(K)
    sort!(σ2,rev=true) # largest to smallest
    σ2 = σ2[1:findlast(σ2.>1e-9)] # get rid of eigenvalues at the end that might be 0 (or below, due to numerical problems)

    talusplot(σ2, isempty(k) ? (1:length(σ2)-1) : k, args...; kwargs...)
end
talusplot(X::Matrix, k::Integer, args...; kwargs...) = talusplot(X, 1:k, args...; kwargs...)


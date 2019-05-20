

# binary search in interval (low,high)
# assumes that if pred(x)==true, then pred(y)==true ∀y>x
# assumes pred(low)==false and pred(high)==true
function binary_search_interval(pred, low::T, high::T, precision::T) where {T<:AbstractFloat}
    while high-low > precision
        mid = (low+high)/2

        if pred(mid)
            high = mid
        else
            low = mid
        end
    end
    (low+high)/2
end
binary_search_interval(pred, low, high, precision) = binary_search_interval(pred, promote(low, high, precision)...)



# Soft thresholding operator S (See Witten et al.)
function soft_threshold(a,c)
    sign(a) .* max( abs(a)-c, 0 )
end





# computes right singular vector of X
function rightsingularvector(X)
    K = Symmetric(X'X);
    N = size(K,1);
    _,v=eig(K,N:N);
    dropdims(v; dims=2)
end



# computes 
#   argmax_u   u'w
#   subject to ||u||_2^2 ≤ 1
# (When used for PCA iterations, w = Xv or w = X'u = (u'X)'.)
function argmax_unpenalized(w)
    # w/norm(w)
    w/max(1e-16,norm(w)) # avoid div by zero
end




# computes 
#   argmax_u   uᵀw
#   subject to ||u||_2^2 ≤ 1
#              uᵀO = (0 0 ... 0)          (vector of zeros)
# (When used for PCA iterations, w = Xv or w = X'u = (u'X)' and U is matrix with previous components, i.e. U and V respectively.)
# See Witten et al 2009, section 3.2.
#     uₖ = Cₖ₋₁Cₖ₋₁ᵀXvₖ / ||Cₖ₋₁ᵀXvₖ||₂
#        = (I - Uₖ₋₁Uₖ₋₁ᵀ)Xvₖ / ||Cₖ₋₁ᵀXvₖ||₂
#        = (Xvₖ - Uₖ₋₁Uₖ₋₁ᵀXvₖ) / √( ||Xvₖ||₂² - ||Uₖ₋₁ᵀXvₖ||₂² )
#        = (w - UUᵀw) / √( ||w||₂² - ||Uᵀw||₂² )
#     where C is the orthogonal complement of O (not needed in computations).
function argmax_unpenalized_orthogonal(w, U)
    # Uw = U'w
    # (w-U*Uw) / sqrt(sum(w.^2)-sum(Uw.^2))

    # Simplified
    u = w - U*(U'w)
    u / norm(u)
end





# computes 
#   argmax_u   u'a
#   subject to ||u||₂² ≤ 1
#              ||u||₁ ≤ c
# (When used for PCA iterations, a = Xv or a = X'u = (u'X)'.)
function argmax_l1penalized(a,c)
    # test with Δ=0 first (i.e. S(a,0)=a)
    # u = a/norm(a)
    u = a/max(1e-16, norm(a)) # avoid div by zero
    c2 = c*c # we almost only use c squared
    (length(a)<=c2 || norm(u,1)≤c) && return u

    # otherwise chose Δ s.t. ||u||₁=c, where u = S(a,Δ)/||S(a,Δ)||₂

    # we could use a heap for b to avoid sorting all of b (should save a lot for very sparse vectors!)
    heap = binary_maxheap(abs(a)) # work with sorted abs(a) since it doesn't change value of Δ and makes the search simpler


    n1 = 0.0
    n2 = 0.0 # squared
    pk=qk = 0.0
    Δ = 0.0
    b = pop!(heap)

    k = 1
    #while k<length(b)
    while k<length(a)
        # pk += b[k]*b[k]
        # qk += b[k]
        pk += b*b
        qk += b

        bNext = pop!(heap)
        # n1Next = qk - k*b[k+1]
        # n2Next = k*b[k+1]*b[k+1] + pk - 2b[k+1]*qk
        n1Next = qk - k*bNext
        n2Next = k*bNext*bNext + pk - 2bNext*qk

        n1Next*n1Next≥c2*n2Next && break # squared version, ok since everything is positive - now we know bₖ₋₁≥Δ≥bₖ        

        n1,n2 = n1Next,n2Next
        b = bNext
        k+=1
    end
    
    A = (k*k-k*c2)
    B = 2(k*n1 - c2*n1)
    C = n1*n1-c2*n2

    p = B/A
    q = C/A

    #Δ = max(0.0, b[k]-(-p/2+sqrt(max(p*p/4-q,0.0)))) # max with zero only to avoid numerical problems
    Δ = max(0.0, b-(-p/2+sqrt(max(p*p/4-q,0.0)))) # max with zero only to avoid numerical problems


    u=soft_threshold(a,Δ); 
    # u/norm(u)
    u/max(1e-16, norm(u)) # avoid div by zero
end




# numerically compute rank 1 approximation
# 2nd component are retrieved using pca_rank(X-u*s*v',...), where (u,s,v) is the result from the first, and so on.
function pcarank1l1(X,c1,c2,tol,maxIters)
    P = size(X,1)
    N = size(X,2)

    v = rightsingularvector(X) # start with optimal v for non-penalized problem.
    u = zeros(P)
    s = 0.0
    converged = false

    # for i=1:maxIters # TODO: convergence criterion
    #     u = c1==Inf ? argmax_unpenalized(X*v) : argmax_l1penalized(X*v,c1)
    #     v = c2==Inf ? argmax_unpenalized(X'u) : argmax_l1penalized(X'u,c2)
    # end
    # s = (u'X*v)[1]


    for i=1:maxIters
        prevU = u
        u = c1==Inf ? argmax_unpenalized(X*v) : argmax_l1penalized(X*v,c1)
        
        Xu = X'u
        prevV = v
        v = c2==Inf ? argmax_unpenalized(Xu) : argmax_l1penalized(Xu,c2)

        prevS = s
        s = (v'Xu)[1]
        #println("s: ", s)

        #println("i: ", i, ", v: ", norm(v-prevV), ", u: ", norm(u-prevU), ", s:", s-prevS)

        norm(v-prevV)<tol && norm(u-prevU)<tol && (converged=true; break)
        s-prevS < tol && (converged=true; break)

        #norm(v-prevV)<tol && norm(u-prevU)<tol && (converged=true; println(norm(v-prevV), " ", norm(u-prevU)); break)
        #s-prevS < tol && (converged=true; println(s-prevS); break)
    end

    (u,s,v,converged)
end






# numerically compute rank 1 approximation
# 2nd component is retrieved using pcarank1l1_orthogonal(X-u*s*v',c1,V), where (u,s,v) is the result from the first, and so on.
function pcarank1l1_orthogonal(X,c1,V,v0,tol,maxIters)
    P = size(X,1)
    N = size(X,2)

    v = v0    

    u = zeros(P)
    s = 0.0
    converged = false

    for i=1:maxIters
        prevU = u
        u = c1==Inf ? argmax_unpenalized(X*v) : argmax_l1penalized(X*v,c1)
        
        Xu = X'u
        prevV = v
        v = isempty(V) ? argmax_unpenalized(Xu) : argmax_unpenalized_orthogonal(Xu,V)

        prevS = s
        s = (v'Xu)[1]

        norm(v-prevV)<tol && norm(u-prevU)<tol && (converged=true; break)
        s-prevS < tol && (converged=true; break)
    end

    (u,s,v,converged)
end


function pcarank1l1_orthogonal(X,c1,V,tol,maxIters)
    # v0 = rightsingularvector(X) # start with optimal v for non-penalized problem. 

    # This enforces orthogonality to V for initial guess.
    v0 = isempty(V) ? rightsingularvector(X) :
                      rightsingularvector( X - (X*V)*V' )

    pcarank1l1_orthogonal(X,c1,V,v0,tol,maxIters)
end





# Sparce PCA with l1-constraints on U-side and orthogonality constraints on V-side.
# If a component with higher optimum than a previous component is found, we backtrack and use it as a starting guess for the previous component.
# This guarantees that S will be sorted and reduces some issues with getting stuck in local optima.
function sparsepcal1(X, k::Integer, c1; backtrack=true, tol=1e-6, maxIters=1000)
    P = size(X,1)
    N = size(X,2)
    @assert N<=P # computations are arranged to exploit N<P (or rather, N<<P)

    
    # set L1-constraints to Inf if they are not active (better for pcarank1l1)
    c1 = c1<√P ? c1 : Inf

    U = zeros(P,k)
    S = zeros(k)
    V = zeros(N,k)

    XOrig = backtrack ? copy(X) : similar(X,0,0) # only necessary to keep track of original matrix if we are backtracking

    v0 = rightsingularvector(X)


    i=1
    while i<=k
        U[:,i], S[i], V[:,i], converged = pcarank1l1_orthogonal(X,c1,V[:,1:i-1],v0,tol,maxIters)
        converged || @warn("sparsepca component $i did not converge")

        if S[i]<0 || any(abs(V[:,i]'V[:,1:i-1]) .> 1e-6) # In this case, we are not finding any more components. I.e. the rank(XOrig)<i.
            U[:,i], S[i], V[:,i] = .0, .0, .0
            break
        end

        if backtrack && i>1 && S[i]-S[i-1]>tol # Backtrack if we found a component with higher S[i] than the previous one.
            v0 = V[:,i] # Use this optimum as a starting point for the search on the previous matrix. The orthogonality of samples (V) guarantees that the starting point will have the same S value, i.e. we will find a better local optima.
            i = something(findfirst(S[i]-S.>tol),0) - 1 # How far do we need to backtrack?
            X = XOrig - U[:,1:i]*(diagm(S[1:i])*V[:,1:i]') # Restore matrix state s.t. only the i first components have been removed.
        else
            X = X - U[:,i]*(S[i]*V[:,i]') # Remove this component to be able to find next component.
            v0 = rightsingularvector( X - (X*V)*V' ) # Starting guess orthogonal to V.
        end
        i += 1
    end

    # change signs for stability, favouring U's that are positive on average
    sgn = (sum(U,1).>=0)*2.0-1.0
    U .*= sgn
    V .*= sgn

    (U,S,V)
end




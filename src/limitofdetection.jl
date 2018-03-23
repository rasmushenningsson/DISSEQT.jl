

# binary search in interval (low,high)
# assumes that if pred(x)=true, then pred(y)=true ∀y>x
# assumes pred(low)=false and pred(high)=true
function binarysearchinterval{T<:AbstractFloat}(pred, low::T, high::T, precision::T)
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
binarysearchinterval(pred, low, high, precision) = binarysearchinterval(pred, promote(low, high, precision)...)


# binary search predicate function ||log₂(f+α)-log₂(r+α)||_(RMS) < threshold
function lodpred(α, f, r, threshold)
    N = length(f)
    norm(log2.(f.+α)-log2.(r.+α)) <= threshold*sqrt(N)
end
function _limitofdetection{T}(freqsF::Vector{T}, freqsR::Vector{T}, threshold, lgLow, lgHigh)
    @assert size(freqsF)==size(freqsR)
    lodpred(10.0^lgLow,freqsF,freqsR,log2(threshold)) && return 10.0^lgLow # early out
    10.0^binarysearchinterval( α->lodpred(10.0^α,freqsF,freqsR,log2(threshold)), lgLow, lgHigh, 1e-3) # logspace search for α
end


# Two 4x4x4xN arrays, one per strand. N=Number of samples.
function _limitofdetection{T}(freqs::Array{T,4}, freqsF::Array{T,4}, freqsR::Array{T,4}, threshold, lgLow, lgHigh)
    N = size(freqs,4)
    @assert size(freqs)==size(freqsF)==size(freqsR)==(4,4,4,N)

    # group samples by consensus codon
    consensusInd = map(i->indmax(freqs[:,:,:,i]), 1:N)

    lod = (10.0^lgLow)*ones(T,4,4,4)

    # compute per group and per variant
    for c in unique(consensusInd)
        mask = consensusInd.==c
        fF = freqsF[:,:,:,mask]
        fR = freqsR[:,:,:,mask]
        for k=1:4,j=1:4,i=1:4
            x = _limitofdetection(fF[i,j,k,:][:], fR[i,j,k,:][:], threshold, lgLow, lgHigh)
            lod[i,j,k] = max(x, lod[i,j,k])
        end
    end
    lod
end




# All swarms should be from the same run.
function _limitofdetection(swarms::AnnotatedArray, swarmsF::AnnotatedArray, swarmsR::AnnotatedArray, threshold, lgLow, lgHigh)
    @assert size(swarms)==size(swarmsF)==size(swarmsR)
    @assert swarms[:id]==swarmsF[:id]==swarmsR[:id]
    P = size(swarms,4)
    N = size(swarms,5)

    lod = zeros(4,4,4,P)
    # TODO: parallelize here???
    for p=1:P
        lod[:,:,:,p] = _limitofdetection(convert(Array,swarms)[:,:,:,p,:],
                                         convert(Array,swarmsF)[:,:,:,p,:],
                                         convert(Array,swarmsR)[:,:,:,p,:],
                                         threshold, lgLow, lgHigh)
    end
    lod
end


function limitofdetection(swarms::AnnotatedArray, swarmsF::AnnotatedArray, swarmsR::AnnotatedArray; groupID=ones(Int,size(swarms,5)), threshold=1.5, lgLow=-6, lgHigh=1)
    @assert length(groupID)==size(swarms,5)

    P = size(swarms,4)
    uniqueGroups = unique(groupID)

    lods = zeros(4,4,4,P,length(uniqueGroups))
    # TODO: parallelize here???
    for (i,g) in enumerate(uniqueGroups)
        mask = groupID.==g
        lods[:,:,:,:,i] = _limitofdetection(swarms[:,:,:,:,mask], swarmsF[:,:,:,:,mask], swarmsR[:,:,:,:,mask], threshold, lgLow, lgHigh)
    end
    lod = squeeze(maximum(lods, 5),5)
    lod = AnnotatedArray(lod)
    annotate!(lod, :segment,   squeeze(swarms[:segment],5))
    annotate!(lod, :position,  squeeze(swarms[:position],5))
    annotate!(lod, :n1,        squeeze(swarms[:n1],5))
    annotate!(lod, :n2,        squeeze(swarms[:n2],5))
    annotate!(lod, :n3,        squeeze(swarms[:n3],5))
    annotate!(lod, :aminoacid, squeeze(swarms[:aminoacid],5))
    annotate!(lod, :codon,     squeeze(swarms[:codon],5))
    lod
end



_featurepositions(a::Integer, b::Integer) = a:3:b
function _featurepositions(a::AbstractString, b::AbstractString)
    a = Int[ parse(Int,x) for x in split(a,',') ]
    b = Int[ parse(Int,x) for x in split(b,',') ]
    vcat( map((x,y)->x:3:y, a, b)... )
end


function featurepositions(featureTable::DataFrame, featureName)
    row = findfirst(isequal(featureName), featureTable[!,:Feature])
    @assert row!=0
    _featurepositions(featureTable[row,:Start], featureTable[row,:End])
end

function featuresegment(featureTable::DataFrame, featureName)
    row = findfirst(isequal(featureName), featureTable[!,:Feature])
    @assert row!=0
    featureTable[row,:Segment]
end




# returns the most specific (i.e. shortest) feature covering the given position
function featureat(featureTable::DataFrame, position::Integer)
    ind = findall(featureTable[!,:Start].<=position.<=featureTable[!,:End]) # matching features
    isempty(ind) && return ""
    i = argmin( featureTable[ind,:End] - featureTable[ind,:Start] )
    featureTable[ind[i],:Feature] # name of shortest feature
end

function codonposition(featureTable::DataFrame, featureName, position::Integer)
    i = findfirst(isequal(featureName), featureTable[!,:Feature])
    @assert i>0 "Feature $featureName not found in table"
    @assert featureTable[i,:Start]<=position<=featureTable[i,:End]
    p = position-featureTable[i,:Start]
    @assert mod(p,3)==0 "position does does not match the ORF"
    div(p,3)+1
end


function annotatevariants!(swarm::AnnotatedArray, featureTable::DataFrame, reference::AbstractString)
    pos = swarm[:position]
    feature = map(p->featureat(featureTable,p), pos)
    refAA = codon2aa(seq2codons(reference))[pos]
    fpos = map((f,p)->codonposition(featureTable,f,p), feature, pos)
    variantAA = swarm[:aminoacid]

    sz = map(max, size(pos), size(variantAA)) # replace with broadcast_shape?
    # variant in style: feature-refAA-featurePos-variantAA
    variant = broadcast!( (f,r,p,v)->"$f-$r$p$v", Vector{String}(undef,sz), feature, refAA, fpos, variantAA)

    annotate!(swarm,:feature,feature)
    annotate!(swarm,:variant,variant)
end



function variantdescriptions(pos, codons, featureTable::DataFrame, reference::AbstractString)
    pos    = reshape(pos,1,length(pos))
    codons = reshape(codons,length(codons),1)

    feature = map(p->featureat(featureTable,p), pos)
    refAA = codon2aa(seq2codons(reference))[pos]
    fpos = map((f,p)->codonposition(featureTable,f,p), feature, pos)
    variantAA = codon2aa(codons)

    sz = map(max, size(pos), size(variantAA)) # replace with broadcast_shape?
    # variant in style: feature-refAA-featurePos-variantAA
    variants = broadcast!( (f,r,p,v)->"$f-$r$p$v", Vector{String}(undef,sz), feature, refAA, fpos, variantAA)
    variants[:]
end

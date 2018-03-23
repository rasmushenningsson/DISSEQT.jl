using BamReader


const NucQualityDict = Dict{UInt8,Int}


# dimensions: nuc, pos
# where nuc is A,C,G or T
const NucQualityMap = Array{NucQualityDict,2}

# init Dict for each nuc and position
nucqualitymap(len::Int) = NucQualityDict[NucQualityDict() for n=1:4, i=1:len]


const NucQualityAccumulator = Vector{NucQualityMap}

nucqualityaccumulator(lengths::Vector{Int}) = NucQualityMap[nucqualitymap(l) for l in lengths]
nucqualityaccumulator(len::Int) = nucqualityaccumulator([len])


# increase the count of a nucleotide quality entry
function inc!(m::NucQualityMap,refPos::Int,nucInd::Int,q::UInt8)
	d = m[nucInd,refPos] # the dict for this position and nucleotide
	d[q] = get!(d,q,0)+1 # is there a better way to increase a value for a key?
end
function inc!(a::NucQualityAccumulator,refInd::Int,refPos::Int,nucInd::Int,q::UInt8)
	inc!(a[refInd],refPos,nucInd,q)
end



# qualityfilter! replaces all quality values below threshold with 0
function qualityfilter!(d::NucQualityDict, threshold::Integer, scratch::Vector{UInt8}=Vector{UInt8}())
	isempty(d) && return
	resize!(scratch,0)

	# get all that should be filtered
	for k in keys(d)
		if k<threshold
			push!(scratch,k)
		end
	end

	# for each filtered entry
	for k in scratch	
		pop!(d,k)
	end
end
function qualityfilter!(m::NucQualityMap, threshold::Integer, scratch::Vector{UInt8}=Vector{UInt8}())
	for d in m
		qualityfilter!(d,threshold,scratch)
	end
end
function qualityfilter!(a::NucQualityAccumulator, threshold::Integer, scratch::Vector{UInt8}=Vector{UInt8}())
	for m in a
		qualityfilter!(m,threshold,scratch)
	end
end





function processread!(nqa::NucQualityAccumulator, r::BamRead)
	refInd = refID(r)+1 # 0 to 1-based index
	
	for c in r
		o = op(c)

		if o==CIGAR_M || o==CIGAR_X || o==CIGAR_Eq # the only chunks that are aligned to both read and reference
			refPos = refpos(c)
			cSeq = seq(r,c)
			cQual = qual(r,c)

			# for all nucleotides in this read chunk
			for (n,q) in zip(cSeq,cQual)
				nuc = base_index(n)
				nuc != 0 && inc!(nqa,refInd,refPos,nuc,q)
				refPos += 1
			end
		end
	end
end





# strands is one of: :both, :forward, :reverse
function nucqualitycount(bamFile::BamFile; strands=:both, mappingQualityThreshold::Int=30, ignoreChimeric::Bool=true)
	@assert strands∈[:both, :forward, :reverse] "Strand must be one of :both, :forward, :reverse."

	seqs = sequences(bamFile)
	seqLengths = Int[s[2] for s in seqs]

	flagsUnset = BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP
	flagsSet = 0

	strands == :forward && (flagsUnset |= BAM_FREVERSE)
	strands == :reverse && (flagsSet   |= BAM_FREVERSE)


	nqa = nucqualityaccumulator(seqLengths)
	processbam!(bamFile,nqa, flagsSet=flagsSet, flagsUnset=flagsUnset,
		        mappingQualityThreshold=mappingQualityThreshold,
		        ignoreChimeric=ignoreChimeric)
	nqa
end




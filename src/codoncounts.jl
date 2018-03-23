# const codons = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"] 


const CodonQualityTriplet = Tuple{UInt8,UInt8,UInt8}
const CodonIndex = Tuple{Int,Int,Int}

codonqualitytriplet() = (0x00,0x00,0x00)::CodonQualityTriplet
codonindex() = (0,0,0)::CodonIndex


valid(ci::CodonIndex) = ci[1]>0 && ci[2]>0 && ci[3]>0
shiftin{T}(t::Tuple{T,T,T},s::T) = (t[2],t[3],s)

function codon(ci::CodonIndex)
	const n = ['A','C','G','T']
	string(ci[1]>0?n[ci[1]]:'x',
		   ci[2]>0?n[ci[2]]:'x',
		   ci[3]>0?n[ci[3]]:'x')
end



const CodonQualityDict = Dict{CodonQualityTriplet,Int}


# dimensions: c, b, a, pos
# where abc is a codon
const CodonQualityMap = Array{CodonQualityDict,4}

# init Dict for each codon and position
codonqualitymap(len::Int) = CodonQualityDict[CodonQualityDict() for n1=1:4,n2=1:4,n3=1:4, i=1:len]


const CodonQualityAccumulator = Array{CodonQualityMap,1}

codonqualityaccumulator(lengths::Array{Int,1}) = CodonQualityMap[codonqualitymap(l) for l in lengths]
codonqualityaccumulator(len::Int) = codonqualityaccumulator([len])


# increase the count of a codon quality triplet
function inc!(m::CodonQualityMap,refPos::Int,codonInd::CodonIndex,
			  q::CodonQualityTriplet)
	d = m[codonInd[3],codonInd[2],codonInd[1],refPos] # the dict for this position and codon
	d[q] = get!(d,q,0)+1 # is there a better way to increase a value for a key?
end
function inc!(a::CodonQualityAccumulator,refInd::Int,refPos::Int,
			  codonInd::CodonIndex,q::CodonQualityTriplet)
	inc!(a[refInd],refPos,codonInd,q)
end



# qualityfilter! replaces all quality values below threshold with 0
#                and merges keys that are equal
function qualityfilter!(d::CodonQualityDict, threshold::Integer, scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
	isempty(d) && return

	resize!(scratch,0)

	# get all that should be filtered
	#for (k,v) in d
	for k in keys(d)
		if k[1]<threshold || k[2]<threshold || k[3]<threshold
			push!(scratch,k)
		end
	end

	# for each filtered entry
	for k in scratch	
		v = pop!(d,k) # remove and get value

		# update quality triplet
		qt = (k[1]>=threshold ? k[1] : 0x00,
			  k[2]>=threshold ? k[2] : 0x00,
			  k[3]>=threshold ? k[3] : 0x00)

		if qt[1]>0 || qt[2]>0 || qt[3]>0
			d[qt] = get!(d,qt,0) + v
		end
	end
end
function qualityfilter!(m::CodonQualityMap, threshold::Integer, scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
	for d in m
		qualityfilter!(d,threshold,scratch)
	end
end
function qualityfilter!(a::CodonQualityAccumulator, threshold::Integer, scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
	for m in a
		qualityfilter!(m,threshold,scratch)
	end
end


iscomplete(qt::CodonQualityTriplet) = qt[1]>0 && qt[2]>0 && qt[3]>0
hascompletecodon(d::CodonQualityDict) = any(iscomplete,keys(d))
function nbrcompletecodons(d::CodonQualityDict)
	s = 0
	for (k,v) in d
		iscomplete(k) && (s+=v)
	end
	s
end


# function removeambiguous!(d::CodonQualityDict,complete::Array{Bool,3},scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
# end

# handles 4x4x4 matrix for a given position in the reference
function removeambiguous!{T}(a::AbstractArray{T,3},scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
	complete = map(hascompletecodon,a)

	for c1=1:4
		for c2=1:4
			for c3=1:4
				d = a[c3,c2,c1]

				resize!(scratch,0)
				for k in keys(d)

					# TODO: code is not type stable (Colon/Int for indexing), fix 
					if !any(complete[ k[3]==0 ? (:) : c3, 
						              k[2]==0 ? (:) : c2, 
						              k[1]==0 ? (:) : c1 ])
						#println((c3,c2,c1))
						push!(scratch,k)
					end
				end

				for k in scratch
					delete!(d,k)
				end
			end
		end
	end
end

function removeambiguous!(m::CodonQualityMap,scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
	for p=1:size(m,4)
		#println(p)
		removeambiguous!(view(m,:,:,:,p),scratch)
	end
end
function removeambiguous!(a::CodonQualityAccumulator,scratch::Array{CodonQualityTriplet,1}=Array{CodonQualityTriplet,1}())
	for m in a
		removeambiguous!(m,scratch)
	end
end




function processread!(cqa::CodonQualityAccumulator, r::BamRead)
	# println(pos(r), " ", cigar(r))
	# println(seq(r))

	refInd = refID(r)+1 # 0 to 1-based index
	
	refPos = -1
	codonIndex = codonindex()
	qt = codonqualitytriplet()

	for c in r
		o = op(c)

		if o==CIGAR_M || o==CIGAR_X || o==CIGAR_Eq # the only chunks that are aligned to both read and reference

			# reset counter if the last chunk didn't end where the current start
			p = refpos(c)
			refPos+1 != p && (codonIndex=codonindex()) # invalidate codon
			refPos = p

			cSeq = seq(r,c)
			cQual = qual(r,c)

			for (i,n,q) in zip(1:len(c),cSeq,cQual)				
				codonIndex = shiftin(codonIndex,base_index(n)) # anything but A=1,C=2,G=3,T=4 will be 0 and invalidate
				qt = shiftin(qt,q)

				if valid(codonIndex)
					# refPos-2 to get first nuc in codon
					inc!(cqa,refInd,refPos-2,codonIndex,qt)
				end

				refPos += 1
			end
		else
			# not a chunk that is aligned to both read and reference
			refPos = -1 # reset pos since we have to start over on a new codon
		end
	end
end




# strands is one of: :both, :forward, :reverse
function codonqualitycount(bamFile::BamFile; strands=:both, mappingQualityThreshold::Int=30, ignoreChimeric::Bool=true)
	@assert strands∈[:both, :forward, :reverse] "Strand must be one of :both, :forward, :reverse."

	seqs = sequences(bamFile)
	seqLengths = Int[s[2] for s in seqs]

	flagsUnset = BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP
	flagsSet = 0

	strands == :forward && (flagsUnset |= BAM_FREVERSE)
	strands == :reverse && (flagsSet   |= BAM_FREVERSE)


	cqa = codonqualityaccumulator(seqLengths)
	processbam!(bamFile,cqa, flagsSet=flagsSet, flagsUnset=flagsUnset,
		        mappingQualityThreshold=mappingQualityThreshold,
		        ignoreChimeric=ignoreChimeric)
	cqa
end



# test function to compare with older basecountextractor results
# TODO: remove
function savetxt(filename::String,cqa::CodonQualityAccumulator,refInd=1)
	os = open(filename,"w")
	m = cqa.list[refInd]

	for refPos=1:size(m.map,4)
		for c1=1:4
			for c2=1:4
				for c3=1:4
					# output dict
					d = m[refPos,(c1,c2,c3)]

					dKeys = collect(keys(d))
					dVals = collect(values(d))

					ind = sortperm(dKeys)
					dKeys = dKeys[ind]
					dVals = dVals[ind]

					for k in dKeys
						print(os, Int(k[1]), ' ', Int(k[2]), ' ', Int(k[3]), ' ')
					end
					println(os)
					for v in dVals
						print(os, v, ' ')
					end
					println(os)


					# wrong order
					# for k in keys(d)
					# 	print(os, Int(k[1]), ' ', Int(k[2]), ' ', Int(k[3]), ' ')
					# end
					# println(os)
					# for v in values(d)
					# 	print(os, v, ' ')
					# end
					# println(os)
				end
			end
		end
	end


	close(os)
end





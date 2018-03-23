using BamReader


type ConsensusInfo
	nucCounts::Array{Float64,2} # counts for A,C,G,T
	nextBaseAligned::Array{Int,1} # counts for consecutive bases by covered by the same read
	deletions::Array{Dict{Int,Int},1} # index i means deletion after i (length->reads supporting)
	insertions::Array{Dict{String,Int},1} # index i means insertion after i (string->reads supporting)
end

function ConsensusInfo(refLength::Int)
	nucCounts = zeros(refLength,4)
	nextBaseAligned = zeros(Int,refLength)
	deletions = [Dict{Int,Int}() for i=1:refLength] # to make sure each element is initalized
	insertions = [Dict{String,Int}() for i=1:refLength] # to make sure each element is initalized
	ConsensusInfo(nucCounts,nextBaseAligned,deletions,insertions)
end


function aggregate_reads(bamFile; mappingQualityThreshold::Int=30, ignoreChimeric::Bool=true)
	refSeqs = sequences(bamFile)
	ci = ConsensusInfo[ConsensusInfo(r[2]) for r in refSeqs]



	for (i,r) in enumerate(bamFile)
	# for r in bamFile
		if flagset(r,BAM_FUNMAP) || flagset(r,BAM_FSECONDARY) || flagset(r,BAM_FDUP)
			continue
		end

		mapq(r)<mappingQualityThreshold && continue
		ignoreChimeric && hastag(r,tag"SA") && continue


		ref = refID(r)+1
		s = seq(r)

		for c in r
			o = op(c)

			# TODO: Handle CIGAR X and =. 
			# 	    NB: nextBaseAligned computation gets more messy
			(o==CIGAR_Eq || o==CIGAR_X) && error("CIGAR X and = are currently not suported")

			if o==CIGAR_M
				# aligned chunk

				readOffs = readpos(c)-1
				refOffs  = refpos(c)-1

				j = 1
				while true
					nuc = base_index(s[readOffs+j])
					1<=nuc<=4 && (ci[ref].nucCounts[refOffs+j,nuc] += 1)
					j>=len(c) && break
					ci[ref].nextBaseAligned[refOffs+j] += 1
					j += 1
				end
			elseif o==CIGAR_I
				# insertion
				readOffs = readpos(c)-1
				
				# TODO: implement better way to get substring of seq
				buf = Array{UInt8}(len(c))
				for j=1:len(c)
					buf[j] = base_char(s[readOffs+j])
				end
				str = String(buf)

				refPos = refpos(c)
				list = ci[ref].insertions[refPos-1] # -1 to index insertions by the ref base before the insertion

				list[str] = get(list,str,0)+1

			elseif o==CIGAR_D
				# deletion
				refPos = refpos(c)
				list = ci[ref].deletions[refPos-1] # -1 to index deletions by the ref base before the deletion
				l = len(c)
				list[l] = get(list,l,0)+1
			end
		end
	end

	ci
end




function dict_max(d)
  maxkey, maxvalue = next(d, start(d))[1]
  for (key, value) in d
    if value > maxvalue
      maxkey = key
      maxvalue = value
    end
  end
  maxkey,maxvalue
end
dict_max(d,defaultKey,defaultValue) = isempty(d) ? (defaultKey,defaultValue) : dict_max(d)




function consensus_single(ci::ConsensusInfo, ref, minSupport, indelMinSupport)
	refLength = size(ci.nucCounts,1)
	@assert length(ref)==refLength # given reference must match BAM header length

	io = IOBuffer()

	i = 1
	while i<=refLength

		nuc = indmax(ci.nucCounts[i,:])
		if ci.nucCounts[i,nuc] < minSupport
			print(io,ref[i])
		else
			const bases = ['A','C','G','T']
			print(io,bases[nuc])
		end


		# handle edge (insertion/deletion/nothing)
		
		maxInsertion, maxInsertionCount = dict_max(ci.insertions[i],"",0)
		maxDeletion, maxDeletionCount   = dict_max(ci.deletions[i],0,0)

		edgeThresh = max(ci.nextBaseAligned[i]+1,indelMinSupport) # >nextBaseAligned, >=indelMinSupport
		if max(maxInsertionCount,maxDeletionCount) >= edgeThresh
			if maxInsertionCount > maxDeletionCount
				print(io,maxInsertion)
			else
				i += maxDeletion # jump over deletion gap
			end
		end

		i += 1
	end
	String(take!(io))
end



function consensus(bamFile::BamFile, fasta::Union{AbstractString,Reference,AbstractArray}=[]; 
				   minSupport::Integer=1, indelMinSupport::Integer=minSupport, 
				   mappingQualityThreshold::Int=30, ignoreChimeric::Bool=true)
	@assert minSupport>0
	@assert indelMinSupport>0
	
	# make sure fasta matches bamFile header
	refSeqs = sequences(bamFile)
	
	if typeof(fasta) <: AbstractString
		fasta = loadfasta(fasta)
	elseif isempty(fasta)
		# create fake fasta with all N's
		fasta = Reference( length(refSeqs) )
		for i=1:length(refSeqs)
			fasta[i] = ( refSeqs[i][1], "N"^refSeqs[i][2] )
		end
	end



	fastaDict = Dict(fasta)
	@assert length(fastaDict)==length(refSeqs) "Reference mismatch between BAM and Fasta"
	for (name,len) in refSeqs
		@assert haskey(fastaDict,name) "Reference sequence \"$name\" not found in fasta"
		@assert length(fastaDict[name])==len "Reference sequence \"$name\" has different length in BAM and Fasta"
	end

	# gather data from BAM
	ci = aggregate_reads(bamFile, mappingQualityThreshold=mappingQualityThreshold, 
		                 ignoreChimeric=ignoreChimeric)


	# make a new fasta file, keeping the same order as in the fasta file
	outFasta = Array{Tuple{String,String},1}(length(fasta))

	for fastaInd = 1:length(fasta)
		# find corresponding index in BAM
		bamInd = findfirst( r -> r[1]==fasta[fastaInd][1], refSeqs )

		c = consensus_single( ci[bamInd], fasta[fastaInd][2], minSupport, indelMinSupport )
		outFasta[fastaInd] = (fasta[fastaInd][1], c)
	end

	outFasta
end

consensus(bamFileName::String, fasta::Union{AbstractString,Reference,AbstractArray}=[]; 
	      minSupport=1, indelMinSupport=minSupport, 
	      mappingQualityThreshold::Int=30, ignoreChimeric::Bool=true) = 
				consensus(BamFile(bamFileName),fasta,minSupport=minSupport,
					      indelMinSupport=indelMinSupport, 
					      mappingQualityThreshold=mappingQualityThreshold, 
					      ignoreChimeric=ignoreChimeric)


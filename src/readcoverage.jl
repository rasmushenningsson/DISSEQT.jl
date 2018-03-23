
type ReadCoverageAccumulator
	dCoverage::Vector{Vector{Int}} # derivative of coverage - makes accumulation faster since only starting/ending points of read chunks must be collected
end

function ReadCoverageAccumulator(seqLengths::Vector{Int})
	ReadCoverageAccumulator(Vector{Int}[zeros(Int,len+1) for len in seqLengths]) # length + 1 to make sure we end at 0
end


function _integratecoverage(x)
	y = cumsum(x)
	@assert y[end]==0
	y[1:end-1]
end

function coverage(rca::ReadCoverageAccumulator)
	N = length(rca.dCoverage)
	map!(_integratecoverage, Vector{Vector{Int}}(N), rca.dCoverage)
end




function processread!(rca::ReadCoverageAccumulator, r::BamRead)
	refInd = refID(r)+1 # 0 to 1-based index
	dCoverage = rca.dCoverage[refInd]
	
	for c in r
		o = op(c)

		if o==CIGAR_M || o==CIGAR_X || o==CIGAR_Eq # the only chunks that are aligned to both read and reference
			pStart = refpos(c)
			dCoverage[pStart] += 1

			pEnd = pStart+len(c)
			dCoverage[pEnd] -= 1
		end
	end
end



function readcoverage(bamFile::BamFile; strands=:both, mappingQualityThreshold::Int=30, ignoreChimeric::Bool=true)
	@assert strands∈[:both, :forward, :reverse] "Strand must be one of :both, :forward, :reverse."

	seqs = sequences(bamFile)
	seqLengths = Int[s[2] for s in seqs]

	flagsUnset = BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP
	flagsSet = 0

	strands == :forward && (flagsUnset |= BAM_FREVERSE)
	strands == :reverse && (flagsSet   |= BAM_FREVERSE)


	rca = ReadCoverageAccumulator(seqLengths)
	processbam!(bamFile, rca, flagsSet=flagsSet, flagsUnset=flagsUnset,
		        mappingQualityThreshold=mappingQualityThreshold,
		        ignoreChimeric=ignoreChimeric)

	coverage(rca)
end

readcoverage(bamFilename::AbstractString; kwargs...) = readcoverage(BamFile(bamFilename);kwargs...)

# get a read coverage matrix per reference sequence
function readcoverage(inputs::AbstractArray; kwargs...) 
	nbrSamples = length(inputs)

	coverage = map(x->readcoverage(x;kwargs...),inputs) # Vector{Vector{Vector{Int}}} : Samples-ReferenceSequence-Coverage

	# reorganize to create one matrix per reference sequence
	nbrSeqsPerSample = map(length, coverage)
	nbrSeqs = nbrSeqsPerSample[1]

	@assert all(nbrSeqsPerSample .== nbrSeqs) "All samples must have the same number of reference sequences."

	covMatrices = Vector{Array{Int,2}}(nbrSeqs)
	for i=1:nbrSeqs
		seqLengthPerSample = map( x->length(x[i]), coverage )
		seqLength = maximum(seqLengthPerSample) # To avoid failing if some samples have indels

		# gather read coverage in matrix filled with zeros at the end
		M = zeros(nbrSamples, seqLength)
		for j=1:nbrSamples
			s = seqLengthPerSample[j]
			M[j,1:s] = coverage[j][i]
		end
		covMatrices[i] = M
	end
	covMatrices
end

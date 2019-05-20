

function reshapesingleton(x::AbstractArray, dim::Integer)
	@assert maximum(size(x))==length(x) "Only one nonsingleton dimension allowed"
	sz = ones(Int,dim)
	sz[dim] = length(x)
	reshape(x,sz...)
end


#function loadswarm(id::AbstractString,freqFolder,consensusFolder="")
function loadswarm(id::AbstractString, freqPath::AbstractString, consensusPath::AbstractString="")
	d = load(freqPath)
	@assert ["codonFreqs", "positions", "coverage"] ⊆ keys(d) "$id.jld does not have the required fields."
	N = length(d["codonFreqs"])

	pos = pop!(d,"positions")
	cf = pop!(d,"codonFreqs")
	cov = pop!(d,"coverage")
	segmentInfo = pop!(d,"segmentInfo")

	# create annotated array - putting segments after each other
	# 4x4x4xP - codons x total number of positions in all segments
	# annotations
		# nucleotide in pos 3, 2, 1
		# codon
		# amino acid
		# position in segment
		# segment
		# coverage
		# consensus
		# annotations from codon frequency computations

	data = cat(cf...; dims=4)


	A = AnnotatedArray(data)

	# id
	# annotate!(A,:id,[id])
	annotate!(A,:id,[id])

	# nucleotides
	nucs = ['A','C','G','T']
	annotate!(A,:n3,copy(nucs))
	annotate!(A,:n2,reshapesingleton(copy(nucs),2))
	annotate!(A,:n1,reshapesingleton(copy(nucs),3))

	# codons
	codons = [ Symbol(String([n1,n2,n3])) for n3 in nucs, n2 in nucs, n1 in nucs ]
	annotate!(A,:codon,codons)

	# amino acids
	aa = codon2aa(codons)
	annotate!(A,:aminoacid,aa)

	
	# nucleotide position
	pos = map( x->reshapesingleton(x,4), pos )
	pos = cat(pos...; dims=4)
	annotate!(A,:position,pos)

	if length(segmentInfo)==1
		annotate!(A,:segment,[Symbol(segmentInfo[1][1])])
	else
		# if multiple segments, we need to annotate each position
		segments = map( x->repmat([Symbol(x[1])], x[2]), segmentInfo )
		segments = vcat(segments...)
		annotate!(A,:segment,reshapesingleton(segments,4))
	end


	# coverage
	cov = map( x->reshapesingleton(x,4), cov )
	cov = cat(cov...; dims=4)
	annotate!(A,:coverage,cov)


	# load consensus if supplied
	if !isempty(consensusPath) 
		# consensus = loadfasta(joinpath(consensusFolder,"$(id)_consensus.fasta"))
		consensus = loadfasta(consensusPath)

		# check that the consensus sequence names and length agrees with the .jld info
		consensusInfo = [(x[1],length(x[2])) for x in consensus]
		@assert segmentInfo==consensusInfo "Mismatch between consensus and codon frequency files"

		# codon consensus
		seqs = map(x->x[2], consensus)
		codonConsensus = map( x->x*"NN", seqs ) # add two NN to the end (we are beyond the reference genome)
		codonConsensus = map(seq2codons, codonConsensus)
		codonConsensus = vcat(codonConsensus...)
		annotate!(A,:consensus,reshapesingleton(codonConsensus,4))
	end



	# other annoations from dict
	for (k,v) in d
		annotate!(A, Symbol(k), [v])
	end

	A
end


function loadswarm(ids::AbstractArray, freqPaths::AbstractArray, consensusPaths::AbstractArray=String[])
	@assert length(ids)==length(freqPaths)
	@assert isempty(consensusPaths) || length(freqPaths)==length(consensusPaths)
	isempty(consensusPaths) && (consensusPaths = repmat([""],length(freqPaths)))
	cat(map((id,f,c)->loadswarm(id,f,c),ids,freqPaths,consensusPaths)...; dims=5)
end


function loadswarm(ids::AbstractArray,freqFolder::AbstractString,consensusFolder::AbstractString="")
	freqPaths = map(id->joinpath(freqFolder,"$id.jld"), ids)
	if isempty(consensusFolder)
		loadswarm(ids, freqPaths)
	else
		loadswarm(ids, freqPaths, map(id->joinpath(consensusFolder,"$(id)_consensus.fasta"), ids))
	end
end






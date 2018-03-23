# for convenience
const Sequence = Tuple{String,String}
const Reference = Array{Sequence,1}


function savefasta(filename::AbstractString, ref::Array{Tuple{String,String},1})
	io = open(filename,"w")

	for (name,seq) in ref
		println(io, '>', name)
		println(io, seq)
	end

	close(io)
end
savefasta(filename::AbstractString, ref::Tuple{String,String}) = savefasta(filename,[ref])


# # FastaIO version
# function loadfasta(ref::AbstractString) 
# 	out = Reference()
# 	for x in FastaReader(ref)
# 		push!(out, (x[1],uppercase(x[2])))
# 	end
# 	out
# end

# Bio.Seq version
function loadfasta(ref::AbstractString) 
	a = collect(open(FASTAReader{DNASequence},ref))
	Sequence[(seqname(s),convert(String,sequence(s))) for s in a]
end





function seq2codons(reference::AbstractString)
    Symbol[reference[i:i+2] for i in 1:length(reference)-2]
end
seq2codons{S<:AbstractString,T<:AbstractString}(reference::Tuple{S,T}) = (reference[1], seq2codons(reference[2]))
seq2codons(references::AbstractArray) = map(seq2codons, references)



const codon2aaDict = Dict{Symbol,Symbol}(:AAA=>:K, :AAC=>:N, :AAG=>:K, :AAT=>:N, :ACA=>:T, :ACC=>:T, :ACG=>:T, :ACT=>:T, :AGA=>:R, :AGC=>:S, :AGG=>:R, :AGT=>:S, :ATA=>:I, :ATC=>:I, :ATG=>:M, :ATT=>:I, :CAA=>:Q, :CAC=>:H, :CAG=>:Q, :CAT=>:H, :CCA=>:P, :CCC=>:P, :CCG=>:P, :CCT=>:P, :CGA=>:R, :CGC=>:R, :CGG=>:R, :CGT=>:R, :CTA=>:L, :CTC=>:L, :CTG=>:L, :CTT=>:L, :GAA=>:E, :GAC=>:D, :GAG=>:E, :GAT=>:D, :GCA=>:A, :GCC=>:A, :GCG=>:A, :GCT=>:A, :GGA=>:G, :GGC=>:G, :GGG=>:G, :GGT=>:G, :GTA=>:V, :GTC=>:V, :GTG=>:V, :GTT=>:V, :TAA=>:*, :TAC=>:Y, :TAG=>:*, :TAT=>:Y, :TCA=>:S, :TCC=>:S, :TCG=>:S, :TCT=>:S, :TGA=>:*, :TGC=>:C, :TGG=>:W, :TGT=>:C, :TTA=>:L, :TTC=>:F, :TTG=>:L, :TTT=>:F)

codon2aa(c::Symbol) = get(codon2aaDict,c,:?)
codon2aa(c::AbstractString) = codon2aa(convert(Symbol,c))

codon2aa(A::AbstractArray) = map!(codon2aa,Array{Symbol}(size(A)),A)

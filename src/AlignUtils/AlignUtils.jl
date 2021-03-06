module AlignUtils

using Distributed
using DataFrames
using Gadfly
using Colors
using JLD
# using Levenshtein
using ..Levenshtein # TODO: use the above line when the Levenshtein package is fixed.
using SynapseClient
using ..BamReader
using ..DISSEQT
using ..SynapseTools


export
	Sample,
	find_samples,
	find_aligned,
	assign_reference!,
	align_sample!,
	align_samples!,
	askforconfirmation,
	makecleanfolder,
	getreferenceinfo,
	reference_sanity_check,
	computecodonfrequencies,
	computenucleotidefrequencies,
	uploadaligned,
	uploadswarms,
	coverageplots

include("log.jl")
include("sample.jl")
include("synapseutils.jl")
include("align.jl")
include("misc.jl")
include("consensus.jl")
include("codonfrequency.jl")
include("nucleotidefrequency.jl")
include("readcoverage.jl")

end

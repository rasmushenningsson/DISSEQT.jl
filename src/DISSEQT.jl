module DISSEQT

# temporary submodule
include("Levenshtein/Levenshtein.jl")

# submodules
include("AnnotatedArrays/AnnotatedArrays.jl")
include("BamReader/BamReader.jl")
include("Kruskal/Kruskal.jl")
include("SynapseTools/SynapseTools.jl")
include("AlignUtils/AlignUtils.jl")
include("Plots/Plots.jl")


using JLD
using StatsBase # for sample
using DataStructures
using DataFrames
using BioSequences

using .AnnotatedArrays
using .BamReader

using .Kruskal


# TODO: get rid of these dependencies - only used by estimatesamplemix
using Convex
using SCS



include("fasta.jl")
include("consensus.jl")
include("bamcounts.jl")
include("codoncounts.jl")
include("codonfrequency.jl")
include("nuccounts.jl")
include("nucfrequency.jl")
include("swarm.jl")
include("timeseries.jl")
include("readcoverage.jl")
include("features.jl")
include("limitofdetection.jl")
include("estimatesamplemix.jl")
include("sparsepca.jl")
include("floydwarshall.jl")
include("mds.jl")
include("hierarchicalclustering.jl")
include("medianfilter.jl")
include("fitnessmodels.jl")
include("isomap.jl")


export 
	Swarm,
	SwarmSegment,
	consensus,
	loadfasta,
	savefasta,
	seq2codons,
	codon2aa,
	codonqualitycount,
	nucqualitycount,
	qualityfilter!,
	removeambiguous!,
	mlcodonfreqs,
	mlnucfreqs,
	positions,
	coverage,
	positionfilter!,
	segments,
	segment,
	segmentinfo,
	segmentnames,
	segmentlengths,
	consensus,
	meta,
	processingmeta,
	loadswarm,
	timeseriesindices,
	timeseries,
	timeseriesat,
	timeseriesmean,
	distancecurve,
	distancecurveequals,
	bifurcationtime,
	evaluatedistancecurve,
	distancecurvemean,
	readcoverage,
	featurepositions,
	featuresegment,
	featureat,
	codonposition,
	annotatevariants!,
	variantdescriptions,
	limitofdetection,
	estimatemix,
	sparsepcal1,
	floydwarshall!,
	floydwarshall,
	mds,
	hierarchicalclustering,
	samplesinclusters,
	medianfilter,
	squareddistances,
	LandscapeModel,
	landscapekernelwidth,
	leaveoneoutlandscapekernelwidth,
	NearestNeighborModel,
	GroupModel,
	predictfitness,
	leaveoneoutpredict,
	varianceunexplained,
	kruskalmds,
	makeconnected!,
	isomapmatrix,
	isomap,
	kruskalisomap




end
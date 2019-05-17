module DISSEQT

using Pkg

haskey(Pkg.installed(),"BamReader") || @warn("Module BamReader not installed. Please refer to DISSEQT installation instructions at https://github.com/rasmushenningsson/DISSEQT.jl")

using JLD
using StatsBase # for sample
using DataStructures
using DataFrames
using BioSequences

using .AnnotatedArrays
using BamReader


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

include("Kruskal.jl")
using .Kruskal

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
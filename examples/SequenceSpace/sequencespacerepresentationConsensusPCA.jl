module SequenceSpaceRepresentationScript

#nprocs()==1 && addprocs() # THREADING

using SynapseClient
using DISSEQT
using DISSEQT.SynapseTools
using DISSEQT.Plots
using DISSEQT.AnnotatedArrays
using LinearAlgebra
using DataFrames
using CSV
using Statistics
using JLD
using Gadfly
using Colors


function cachemetadata(syn,alignmentFolder,metadataID)
    metadataCacheFilename = "metadatacache.csv"
    metadata = DataFrame()
    if isfile(metadataCacheFilename)
        metadata = CSV.read(metadataCacheFilename; missingstrings=["","NA"])
    else
        # Get metadata (metadataID can be the synapse ID of a .csv file or of a folder containing .csv files).
        metadata = getsamplemetadata(syn, metadataID, :Rejected=>"No") # Always remove rejected samples

        # extend metadata
        appendswarmids!(syn,metadata,alignmentFolder)
        appendalignedids!(syn,metadata,alignmentFolder)

        # download swarms and consensus sequences
        downloadswarms!(syn,metadata,cacheDir="MutantSwarms")
        downloadconsensuses!(syn,metadata,cacheDir="MutantSwarms")

        CSV.write(metadataCacheFilename, metadata; missingstring="NA")
    end
    metadata
end


function main()
    syn = SynapseClient.login()

    # projectFolder should point to MyProject folder (containing subfolders Metadata and Analysis)
    projectFolder   = "syn11639899" # VignuzziLabPublic/Projects/FitnessLandscapes
    analysisFolder  = getchildbyname(syn, projectFolder, "Analysis")
    alignmentFolder = getchildbyname(syn, analysisFolder, "Alignment")
    referenceFolder = getchildbyname(syn, alignmentFolder, "ReferenceGenomes")

    isdir("MutantSwarms") || mkdir("MutantSwarms")
    isdir("plots") || mkdir("plots")

    # set to true to upload files
    doUpload = true
    uploadName = "fitness landscapes invitro"

    metadataID = getchildbyname(syn, projectFolder, "Metadata", "fitness landscapes invitro.csv")
    metadata = cachemetadata(syn, alignmentFolder, metadataID)
    metadata = metadata[ .~ismissing.(metadata[:Dose]), : ] # skip samples with unknown dose


    # Get Reference Genomes
    referenceNames = ["WT", "Less", "More", "Stop"]
    referenceIDs, referenceGenomes = getreferencegenomes(syn, referenceFolder, referenceNames)


    # get CVB3 feature list
    featureFileID = getchildbyname(syn,referenceFolder,"CVB3_features.csv")
    features = CSV.read(localpath(syn, featureFileID); missingstrings=["","NA"])
    positions = featurepositions(features, "ORF") # coding positions

    # load swarms
    swarms  = loadswarm(metadata[:SampleID],metadata[:SwarmPath],metadata[:ConsensusPath])

    # filter by positions
    swarms = swarms[:,:,:,inset(:position,positions),:]


    coverage = dropdims(swarms[:coverage]; dims=(1,2,3))
    meanCoverage = mean(coverage; dims=1)[:]

    # filter by mean read coverage at coding sites
    sampleMask = meanCoverage .>= 1000
    println("Removed ", count(.~sampleMask), " samples due to low coverage.")
    swarms  = swarms[:,:,:,:,sampleMask]
    metadata = metadata[sampleMask,:]


    # remove mixed samples
    mixes = estimatemix(swarms, referenceGenomes)
    mixMask = maximum(mixes; dims=2)[:] .>= 0.98
    println("Removed ", count(.~mixMask), " mixed samples.")
    swarms  = swarms[:,:,:,:,mixMask]
    metadata = metadata[mixMask,:]


    println("Total number of samples used: ", size(swarms,5))



    nbrDims = 13 # same as for the SMSSVD representation
    

    # setup dependencies
    dependencies = unique(vcat(featureFileID, referenceIDs, metadata[:SwarmID], metadata[:MetadataID]))


    # get rid of annotations to simplify math stuff
    # arrays are 4x4x4xPxN - nuc3 x nuc2 x nuc1 x positions x samples
    s = convert(Array, swarms)
    s = reshape(s, 64, size(swarms,4), size(swarms,5)) # reshape to VxN - variables x samples
    

    # Convert to consensus (1-of-K) encoding
    M = maximum(s; dims=1)
    M[M.<0.25] .= 1. # if no variant is above 0.25, set the maximum to 1 which will force all variables to zero.
    X = float(s.==M) # 1 if equals to the maximum and 0 otherwise


    # reshape to VxN - variables x samples
    X = reshape(X, 64*size(X,2), size(X,3))

    # center variables
    X = X .- mean(X; dims=2) # remove mean over samples
    


    # PCA keeping only nbrDims dimensions
    N = size(X,2)
    F = eigen(Symmetric(X'X), N-nbrDims+1:N)
    Σ = sqrt.(max.(0.,reverse(F.values)))
    V = F.vectors[:,end:-1:1]
    U = X*V ./ Σ'


    uploads = UploadList()



    # Low-Rank Representation, Sample names, positions and codon list
    seqSpaceRep = "plots/seqspacerepresentationConsensusPCA.jld"
    save(seqSpaceRep, "U",U, "Σ",Σ, "V",V,
                      "SampleID", convert(Vector{String},metadata[:SampleID]),
                      "positions", swarms[:position][:],
                      "codons",    swarms[:codon][:],
                      "consensus", dropdims(swarms[:consensus]; dims=(1,2,3)) )
    markforupload!(uploads, "SequenceSpaceRepresentation", seqSpaceRep, dependencies)


    # make PCA plots

    virusColors = Dict("WT"  =>colorant"black",
                       "Less"=>colorant"blue",
                       "More"=>colorant"green",
                       "Stop"=>colorant"red")
    virusColorMap = groupcolors(metadata[:Reference],virusColors)

    mutagenColors = Dict("5FU" =>RGB(1,0,0),
                         "AML" =>RGB(0,1,0),
                         "AZC" =>RGB(0,0,1),
                         "Mn"  =>RGB(0.8,0.8,0),
                         "Mock"=>RGB(1,0,1),
                         "RIBA"=>RGB(0,1,1))
    mutagenColorMap = groupcolors(metadata[:Mutagen],mutagenColors)

    P = pairwisescatterplot(V, metadata[:Reference], virusColorMap,
                               metadata[:Mutagen],   mutagenColorMap,
                               point_size=0.6mm)
    # markforupload!(uploads, "SequenceSpacePlots", saveplot(P, [:png,:svg,:pdf], "plots/pairwiseConsensusPCA"), dependencies)
    markforupload!(uploads, "SequenceSpacePlots", saveplot(P, [:png,:pdf], "plots/pairwiseConsensusPCA"), dependencies)



    doUpload && uploadfiles(syn, createchildfolder(syn, analysisFolder, "SequenceSpace", uploadName), uploads, @__FILE__, "Sequence Space Representation")


    # swarms, metadata, X, U, Σ, V
    nothing
end

end
SequenceSpaceRepresentationScript.main();



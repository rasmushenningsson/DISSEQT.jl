module TalusPlotScript

using SynapseClient
using DISSEQT
using DISSEQT.Plots
using DISSEQT.AnnotatedArrays
using DISSEQT.SynapseTools
using Statistics
using DataFrames
using CSV
using JLD
using Gadfly
using Cairo
using Fontconfig


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
    mixMask = maximum(mixes;dims=2)[:] .>= 0.98
    println("Removed ", count(.~mixMask), " mixed samples.")
    swarms  = swarms[:,:,:,:,mixMask]
    metadata = metadata[mixMask,:]



    println("Total number of samples used: ", size(swarms,5))


    # load limit of detection
    limitOfDetectionID = getchildbyname(syn,analysisFolder,"LimitOfDetection","fitness landscapes invitro","limitofdetection.jld")
    limitOfDetection = convert(AnnotatedArray,load(localpath(syn,limitOfDetectionID,downloadLocation="MutantSwarms")))
    limitOfDetection = limitOfDetection[:,:,:,inset(:position,positions)] # filter by positions
    @assert isequal(swarms[:position][:], limitOfDetection[:position][:]) # the positions must be the same


    # get rid of annotations to simplify math stuff
    # arrays are 4x4x4xPxN - nuc3 x nuc2 x nuc1 x positions x samples
    s = convert(Array, swarms)
    α = convert(Array, limitOfDetection)

    # Put a lower limit on the limit of detection
    α = max.(α,1e-3)

    # apply transformation
    X = log2.(s.+α)

    # center variables
    X = X .- mean(X;dims=5) # remove mean over samples


    # reshape to VxN - variables x samples
    X = reshape(X, 64*size(X,4), size(X,5))



    # Non-intrusive prefiltering of variables that will never contribute
    σThreshold = 1e-2 # same as the lower limit we use for projection score
    σ = dropdims(std(X;dims=2);dims=2)
    σ = σ/maximum(σ) # normalize
    X = X[σ.>=σThreshold,:]




    uploads = UploadList()
    dependencies = unique(vcat(featureFileID, referenceIDs, limitOfDetectionID, metadata[:SwarmID], metadata[:MetadataID]))



    nbrDims = 40

    pl = talusplot(X, 1:nbrDims, Guide.title("Talus plot"))
    markforupload!(uploads, "plots", saveplot(pl, [:png,:svg,:pdf], "plots/talus", width=29.7cm/2, height=21cm/2), dependencies)

    doUpload && uploadfiles(syn, createchildfolder(syn, analysisFolder, "Dimensionality", uploadName), uploads, @__FILE__, "Talus plot")

    nothing
end

end
TalusPlotScript.main();


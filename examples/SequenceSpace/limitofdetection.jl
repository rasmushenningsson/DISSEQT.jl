module LimitOfDetectionScript

using SynapseClient
using DISSEQT
using DISSEQT.SynapseTools
using DISSEQT.AnnotatedArrays
using Statistics
using DataFrames
using CSV
using JLD


function cachemetadata(syn,alignmentFolder,metadataID)
    metadataCacheFilename = "metadatacache.csv"
    metadata = DataFrame()
    if isfile(metadataCacheFilename)
        metadata = CSV.read(metadataCacheFilename; missingstrings=["","NA"])
    else
        # Get metadata (metadataID can be the synapse ID of a .csv file or of a folder containing .csv files).
        metadata = getsamplemetadata(syn, metadataID, :Rejected=>"No") # Always remove rejected samples

        # extend metadata
        appendswarmids!(syn,metadata,alignmentFolder,[:both,:forward,:reverse])
        appendalignedids!(syn,metadata,alignmentFolder)

        # download swarms and consensus sequences
        downloadswarms!(syn,metadata,cacheDir="MutantSwarms",[:both,:forward,:reverse])
        downloadconsensuses!(syn,metadata,cacheDir="MutantSwarms")

        CSV.write(metadataCacheFilename, metadata; missingstring="NA")
    end
    metadata
end


function main()
    syn = SynapseClient.login()

    # prohectFolder should point to MyProject folder (containing subfolders Metadata and Analysis)
    projectFolder   = "syn11639899" # VignuzziLabPublic/Projects/FitnessLandscapes
    analysisFolder  = getchildbyname(syn, projectFolder, "Analysis")
    alignmentFolder = getchildbyname(syn, analysisFolder, "Alignment")
    referenceFolder = getchildbyname(syn, alignmentFolder, "ReferenceGenomes")

    isdir("MutantSwarms") || mkdir("MutantSwarms")
    isdir("MutantSwarms/forward") || mkdir("MutantSwarms/forward")
    isdir("MutantSwarms/reverse") || mkdir("MutantSwarms/reverse")

    # set to true to upload files
    doUpload = true
    uploadName = "fitness landscapes invitro"

    metadataID = getchildbyname(syn, projectFolder, "Metadata", "fitness landscapes invitro.csv")
    metadata = cachemetadata(syn, alignmentFolder, metadataID)

    # get CVB3 feature list
    featureFileID = getchildbyname(syn,referenceFolder,"CVB3_features.csv")
    features = CSV.read(localpath(syn, featureFileID); missingstrings=["","NA"])

    positions = featurepositions(features, "ORF") # coding positions

    # load swarms
    swarms  = loadswarm(metadata[:SampleID],metadata[:SwarmPath],metadata[:ConsensusPath])
    swarmsF = loadswarm(metadata[:SampleID],metadata[:SwarmForwardPath],metadata[:ConsensusPath])
    swarmsR = loadswarm(metadata[:SampleID],metadata[:SwarmReversePath],metadata[:ConsensusPath])


    coverage = dropdims(swarms[:coverage]; dims=(1,2,3))
    meanCoverage = mean(coverage[positions,:]; dims=1)[:]

    # filter by mean read coverage at coding sites
    sampleMask = meanCoverage .>= 1000
    println("Removed ", count(.~sampleMask), " samples due to low coverage.")
    swarms  = swarms[:,:,:,:,sampleMask]
    swarmsF = swarmsF[:,:,:,:,sampleMask]
    swarmsR = swarmsR[:,:,:,:,sampleMask]
    metadata = metadata[sampleMask,:]


    println("Total number of samples used: ", count(sampleMask))

    
    
    uploads = UploadList()
    dependencies = unique(vcat(featureFileID, metadata[:SwarmID], metadata[:SwarmForwardID], metadata[:SwarmReverseID], metadata[:MetadataID]))


    lod = limitofdetection(swarms, swarmsF, swarmsR, groupID=metadata[:Run])
    lodDict = Dict( string(p[1])=>p[2] for p in pairs(convert(Dict,lod)) )
    lodfn = "limitofdetection.jld"
    save(lodfn, lodDict)
    markforupload!(uploads, "", lodfn, dependencies)


    doUpload && uploadfiles(syn, createchildfolder(syn, analysisFolder, "LimitOfDetection", uploadName), uploads, @__FILE__, "Limit of Detection")


    #swarms,swarmsF,swarmsR,lod,metadata    
    nothing
end

end
LimitOfDetectionScript.main();


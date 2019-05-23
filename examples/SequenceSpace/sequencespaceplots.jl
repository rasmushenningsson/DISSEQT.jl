module SequenceSpacePlotsScript

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



function main()
    syn = SynapseClient.login()

    # projectFolder should point to MyProject folder (containing subfolders Metadata and Analysis)
    projectFolder  = "syn11639899" # VignuzziLabPublic/Projects/FitnessLandscapes
    analysisFolder = getchildbyname(syn, projectFolder, "Analysis")
    fitnessFolder  = getchildbyname(syn, projectFolder, "Fitness")

    isdir("MutantSwarms") || mkdir("MutantSwarms")
    isdir("plots") || mkdir("plots")

    # set to true to upload files
    doUpload = true
    uploadName = "fitness landscapes invitro"

    metadataID = getchildbyname(syn, projectFolder, "Metadata", "fitness landscapes invitro.csv")
    metadata = getsamplemetadata(syn, metadataID, :Rejected=>"No") # Always remove rejected samples
    metadata = metadata[ .~ismissing.(metadata[:Dose]), : ] # skip samples with unknown dose



    seqSpaceRepID = getchildbyname(syn, analysisFolder, "SequenceSpace", uploadName, "SequenceSpaceRepresentation", "seqspacerepresentation.jld")
    seqSpaceRepDict = load(localpath(syn, seqSpaceRepID))
    U = seqSpaceRepDict["U"]
    Σ = seqSpaceRepDict["Σ"]
    V = seqSpaceRepDict["V"]
    sampleID  = seqSpaceRepDict["SampleID"]
    positions = seqSpaceRepDict["positions"]
    codons    = seqSpaceRepDict["codons"]
    consensus = permutedims(seqSpaceRepDict["consensus"],(2,1)) # transpose to get Samples x Codons
    nbrDims = length(Σ)

    # filter metadata by samples included in the Sequence Space Representation
    ind = indexin(sampleID, metadata[:SampleID])
    @assert all(ind.!=0) "Unknown SampleID"
    metadata = metadata[ind,:]
   
    
    Y = V.*Σ' # sample representation scaled by singular values





    # get fitness
    fitnessID = getchildbyname(syn,fitnessFolder,"fitness landscapes invitro fitness.csv")
    appendfitness!(metadata, getsamplefitness(syn, fitnessID))


    # setup dependencies
    dependencies = unique(vcat(seqSpaceRepID, metadata[:MetadataID]))
    fitnessDependencies = unique(vcat(dependencies, collect(skipmissing(metadata[:FitnessID]))))



    uploads = UploadList()



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


    # Isomap plots

    K = Y*Y'
    D = sqrt.(diag(K) .+ diag(K)' .- 2*K) # sample distance matrix
    Z,stress,converged,dists = kruskalisomap(D, 2, 2, maximum(D)/20)
    converged || error("Kruskal's stress MDS did not converge.")

    isomapPlotArgs = (Guide.xlabel("Dim 1"), Guide.ylabel("Dim 2"),
                      Coord.Cartesian(xmin=minimum(Z[:,1]),xmax=maximum(Z[:,1]),ymin=minimum(Z[:,2]),ymax=maximum(Z[:,2])),
                      Theme(key_title_font_size=2*11pt, key_label_font_size=2*8pt, major_label_font_size=2*11pt, minor_label_font_size=2*8pt))

    dfIsomap = DataFrame(x=Z[:,1], y=Z[:,2], Virus=metadata[:Reference], Mutagen=metadata[:Mutagen])
    pl = plot(dfIsomap, x=:x, y=:y, color=:Virus, virusColorMap, isomapPlotArgs...)
    markforupload!(uploads, "SequenceSpacePlots", saveplot(pl, [:png,:svg,:pdf], "plots/isomap_virus"), dependencies)
    pl = plot(dfIsomap, x=:x, y=:y, color=:Mutagen, mutagenColorMap, isomapPlotArgs...)
    markforupload!(uploads, "SequenceSpacePlots", saveplot(pl, [:png,:svg,:pdf], "plots/isomap_mutagen"), dependencies)





    # fitness landscape plots
    fitnessMask = .~ismissing.(metadata[:Fitness])


    fitnessTrain = convert(Vector,metadata[fitnessMask,:Fitness]);
    ZTrain = Z[fitnessMask,:];
    σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(ZTrain.^2;dims=2))) # scale by point size cloud
    @time landscapeModel = LandscapeModel(ZTrain, fitnessTrain, σValues, nbrIter=1000)
    println("Isomap Landscape σ: ", landscapeModel.σ)


    landscapeRes = 500
    xLim = [extrema(Z[:,1])...]
    yLim = [extrema(Z[:,2])...]
    xLim += 0.1*(xLim[2]-xLim[1])*[-1,1] # add margins
    yLim += 0.1*(yLim[2]-yLim[1])*[-1,1] # add margins
    x = range(xLim[1], stop=xLim[2], length=landscapeRes)
    y = range(yLim[1], stop=yLim[2], length=landscapeRes)
    pl = fitnesslandscapeplot(landscapeModel,x,y,Z,metadata[:Reference],virusColors,
                              width=1024, height=768, zAspect=0.15,
                              cameraCenter=(0,0,0), cameraEye=(-.75,.6,.75), markerSize=4,
                              σTransparency=0.05);
    markforupload!(uploads, "FitnessLandscapePlots", saveplot(pl, [:html,:png,:pdf], "plots/isomap_landscape"), fitnessDependencies)
    # saveplot(pl, :html, "plots/isomap_landscape_pyembed", js=:embed) # Save local file with plotlyjs embedded.



    # landscape plot with mutagen colors
    pl = fitnesslandscapeplot(landscapeModel,x,y,Z,metadata[:Mutagen],mutagenColors,
                              width=1024, height=768, zAspect=0.15,
                              cameraCenter=(0,0,0), cameraEye=(-.75,.6,.75), markerSize=4,
                              σTransparency=0.05);
    markforupload!(uploads, "FitnessLandscapePlots", saveplot(pl, [:html,:png,:pdf], "plots/isomap_landscape_mutagen"), fitnessDependencies)
    # saveplot(pl, :html, "plots/isomap_landscape_mutagen_pyembed", js=:embed) # Save local file with plotlyjs embedded.



    doUpload && uploadfiles(syn, createchildfolder(syn, analysisFolder, "SequenceSpace", uploadName), uploads, @__FILE__, "Sequence Space Plots")


    # seqSpaceRepDict, metadata, U, Σ, V, Y, Z
    nothing
end

end
SequenceSpacePlotsScript.main();
# seqSpaceRepDict, metadata, U, Σ, V, Y, Z = SequenceSpacePlotsScript.main();


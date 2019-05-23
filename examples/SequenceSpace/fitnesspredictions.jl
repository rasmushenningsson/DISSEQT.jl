module FitnessPredictionsScript

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



# put mutagen doses into low,medium,high categories
function addlevel!(metadata)
    levelDict = Dict(("5FU","50uM")=>"Low", ("5FU","100uM")=>"Medium", ("5FU","150uM")=>"Medium", ("5FU","200uM")=>"High", ("AML","50uM")=>"Low", ("AML","100uM")=>"Medium", ("AML","200uM")=>"High", ("AZC","50uM")=>"Low", ("AZC","100uM")=>"Medium", ("AZC","200uM")=>"Medium", ("AZC","300uM")=>"High", ("Mn","0.25mM")=>"Medium", ("Mn","0.33mM")=>"Medium", ("Mn","0.5mM")=>"High", ("Mock","0uM")=>"Low", ("RIBA","50uM")=>"Low", ("RIBA","100uM")=>"Medium", ("RIBA","200uM")=>"Medium", ("RIBA","300uM")=>"High")
    metadata[:Level] = map( (m,d)->get(levelDict, (m,d), missing), metadata[:Mutagen], metadata[:Dose] )
end



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
   
    Y = V.*Σ' # sample representation scaled by singular values



    # also load the Consensus-based PCA representation
    seqSpaceRepConsensusID = getchildbyname(syn, analysisFolder, "SequenceSpace", uploadName, "SequenceSpaceRepresentation", "seqspacerepresentationConsensusPCA.jld")
    seqSpaceRepConsensusDict = load(localpath(syn, seqSpaceRepConsensusID))
    UConsensus = seqSpaceRepConsensusDict["U"]
    ΣConsensus = seqSpaceRepConsensusDict["Σ"]
    VConsensus = seqSpaceRepConsensusDict["V"]
    @assert sampleID == seqSpaceRepConsensusDict["SampleID"] # make sure the same samples were used
    nbrDimsConsensus = length(ΣConsensus)

    YConsensus = VConsensus.*ΣConsensus' # sample representation scaled by singular values



    # filter metadata by samples included in the Sequence Space Representation
    ind = indexin(sampleID, metadata[:SampleID])
    @assert all(ind.!=0) "Unknown SampleID"
    metadata = metadata[ind,:]




    addlevel!(metadata) # categorize dosages as low, medium, high


    # get fitness
    fitnessID = getchildbyname(syn,fitnessFolder,"fitness landscapes invitro fitness.csv")
    appendfitness!(metadata, getsamplefitness(syn, fitnessID))






    # setup dependencies
    dependencies = unique(vcat(seqSpaceRepID, seqSpaceRepConsensusID, metadata[:MetadataID], collect(skipmissing(metadata[:FitnessID]))))

    


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


    # Isomap representation
    K = Y*Y'
    D = sqrt.(diag(K) .+ diag(K)' .- 2*K) # sample distance matrix
    Z,stress,converged,dists = kruskalisomap(D, 2, 2, maximum(D)/20)
    converged || error("Kruskal's stress MDS did not converge.")

    # Isomap representation (Consensus PCA)
    KConsensus = YConsensus*YConsensus'
    DConsensus = sqrt.(max.(0., diag(KConsensus) .+ diag(KConsensus)' .- 2*KConsensus)) # sample distance matrix (make sure it's nonnegative to handle rounding errors for colocated points)
    ZConsensus,_,converged,_ = kruskalisomap(DConsensus, 2, 2, maximum(DConsensus)/20)
    converged || error("Kruskal's stress MDS did not converge for Consensus PCA.")




    # prediction
    fitnessMask = .~ismissing.(metadata[:Fitness])


    # make a table with all predictors
    predictionResults = DataFrame(Name=String[], Type=String[], VarianceExplained=Float64[], NbrSamples=Int[])


    # run twice, the second time excluding samples that constitute singleton groups in the "Lineage, Mutagen, Dose" group predictor
    mask = copy(fitnessMask)
    for i=1:2
        fitness = collect(skipmissing(metadata[mask,:Fitness]))

        # Isomap landscape
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(Z[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(Z[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, Z[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] Isomap (2d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )
        # println(predictionResults)

        # Reduced space landscape
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(Y[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(Y[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, Y[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] SMSSVD ($(nbrDims)d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )
        # println(predictionResults)


        # Isomap landscape (Consensus PCA)
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(ZConsensus[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(ZConsensus[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, ZConsensus[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] Consensus Isomap (2d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )
        # println(predictionResults)

        # Reduced space landscape  (Consensus PCA)
        σValues = 10.0.^range(-3,stop=1,length=1000) * sqrt(mean(sum(YConsensus[mask,:].^2; dims=2))) # scale by point size cloud
        @time σPerSample,_,_ = leaveoneoutlandscapekernelwidth(YConsensus[mask,:], fitness, σValues, nbrIter=1000)
        f = leaveoneoutpredict(LandscapeModel, YConsensus[mask,:], fitness, perSampleModelArgs=σPerSample)
        push!(predictionResults, ("[GKS] Consensus PCA ($(nbrDimsConsensus)d)", "Gaussian Kernel Smoother", 1.0.-varianceunexplained(f,fitness), count(mask)) )


        # Nearest neighbor predictor: Reduced dimension
        f = leaveoneoutpredict(NearestNeighborModel, Y[mask,:], fitness)
        push!(predictionResults, ("[NN] SMSSVD ($(nbrDims)d)", "Nearest Neighbor", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Nearest neighbor predictor: Consensus
        f = leaveoneoutpredict(NearestNeighborModel, consensus[mask,:], fitness)
        push!(predictionResults, ("[NN] Consensus", "Nearest Neighbor", 1.0.-varianceunexplained(f,fitness), count(mask)) )


        # Group predictor: Lineage, Drug
        f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Mutagen]], fitness)
        push!(predictionResults, ("[G] Lineage/Mutagen", "Group", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Group predictor: Lineage, Level
        f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Level]], fitness)
        push!(predictionResults, ("[G] Lineage/Dose", "Group", 1.0.-varianceunexplained(f,fitness), count(mask)) )

        # Group predictor: Lineage, Drug, Dose
        #f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Mutagen,:Dose]], fitness)
        f = leaveoneoutpredict(GroupModel, metadata[mask,[:Reference,:Mutagen,:Level]], fitness)
        push!(predictionResults, ("[G] Lineage/Mutagen/Dose", "Group", 1.0.-varianceunexplained(f,fitness), count(mask)) )
        newMask = i==1 ? .~isnan.(f) : trues(size(f)) # For next iteration, exclude singleton groups
        i==1 && println(count(.~newMask), " singleton groups removed.")

        mask[mask] = newMask
    end

    println(predictionResults)
    CSV.write("plots/fitnesspredictions.csv", predictionResults; missingstring="NA")
    markforupload!(uploads, "FitnessPredictionPlots", "plots/fitnesspredictions.csv", dependencies)

    N = div(size(predictionResults,1),2)
    xMax = maximum(predictionResults[.~isnan.(predictionResults[:VarianceExplained]), :VarianceExplained])

    # with all samples
    df = predictionResults[1:N,:]
    df = df[.~isnan.(df[:VarianceExplained]),:] # remove predictors that couldn't be used ([G] Lineage/Mutagen/Dose)
    plotArgs = (layer(df, x=:VarianceExplained, color=:Name, Geom.bar(orientation=:horizontal)),
                Coord.Cartesian(xmin=0, xmax=xMax, yflip=true),
                Guide.yticks(ticks=nothing),
                Guide.xlabel("Variance Explained"), Guide.ylabel(""))
    themeArgs = (:key_title_font_size=>2*11pt, :key_label_font_size=>2*8pt, :major_label_font_size=>2*11pt, :minor_label_font_size=>2*8pt, :point_label_font_size=>2*8pt)

    pl = plot(plotArgs..., Theme(;themeArgs...))
    markforupload!(uploads, "FitnessPredictionPlots", saveplot(pl, [:png,:svg,:pdf], "plots/fitnesspredictions_allsamples"), dependencies)

    pl = plot(plotArgs..., Theme(key_position=:none; themeArgs...)) # no legend
    markforupload!(uploads, "FitnessPredictionPlots", saveplot(pl, [:png,:svg,:pdf], "plots/fitnesspredictions_allsamples_nolegend"), dependencies)



    # with singleton groups removed
    df = predictionResults[N+1:end,:]
    df = df[.~isnan.(df[:VarianceExplained]),:] # remove predictors that couldn't be used (none)
    plotArgs = (layer(df, x=:VarianceExplained, color=:Name, Geom.bar(orientation=:horizontal)),
                Coord.Cartesian(xmin=0, xmax=xMax, yflip=true),
                Guide.yticks(ticks=nothing),
                Guide.xlabel("Variance Explained"), Guide.ylabel(""))
    themeArgs = (:key_title_font_size=>2*11pt, :key_label_font_size=>2*8pt, :major_label_font_size=>2*11pt, :minor_label_font_size=>2*8pt, :point_label_font_size=>2*8pt)

    pl = plot(plotArgs..., Theme(;themeArgs...))
    markforupload!(uploads, "FitnessPredictionPlots", saveplot(pl, [:png,:svg,:pdf], "plots/fitnesspredictions"), dependencies)

    pl = plot(plotArgs..., Theme(key_position=:none; themeArgs...)) # no legend
    markforupload!(uploads, "FitnessPredictionPlots", saveplot(pl, [:png,:svg,:pdf], "plots/fitnesspredictions_nolegend"), dependencies)



    
    doUpload && uploadfiles(syn, createchildfolder(syn, analysisFolder, "SequenceSpace", uploadName), uploads, @__FILE__, "Sequence Space Plots")


    # seqSpaceRepDict, metadata, U, Σ, V, Y, Z
    nothing
end

end
FitnessPredictionsScript.main();
# seqSpaceRepDict, metadata, U, Σ, V, Y, Z = FitnessPredictionsScript.main();


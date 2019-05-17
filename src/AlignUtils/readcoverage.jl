
import SynapseClient: Activity


function _savecoverageplot(layer,format,outPath,args...)
    if format==:png
        args = (Theme(background_color=colorant"white"), args...)
    end

    pl = plot(layer,args...)

    filename = string(outPath,'.',format)
    if format==:png
        draw(PNG(filename, 29.7cm, 21cm), pl)
    elseif format==:pdf
        draw(PDF(filename, 29.7cm, 21cm), pl)
    elseif format==:svg
        draw(SVG(filename, 29.7cm, 21cm), pl)
    else
        error("Unknown image format $format")
    end
    filename
end


# outputs list of files created
function _coverageplots(samplePaths, sampleNames, outputFolder, outName, outFormats, log)
    print(log, "Computing read coverage...")
    coverage = readcoverage(samplePaths)
    println(log, "Done.")

    # check if these samples have one or multiple segments
    nbrSegments = length(coverage)

    segmentNames = [""]
    if nbrSegments>1
        seqs = sequences(BamFile(samplePaths[1])) # array of (segmentName,length)
        segmentNames = String[x[1] for x in seqs]
    end

    outFiles = Vector{String}()
    outSegment = Vector{String}()

    print(log, "Plotting...")
    for (i,segmentName) in enumerate(segmentNames)
        # Create DataFrame with all the data for the plot
        cov = coverage[i]
        maximum(cov)==0 && (cov[1]=1) # Gadfly fix for logscale when no sample has data.
        nbrSamples = size(cov,1)
        nbrPos     = size(cov,2)

        pos = repmat( (1:nbrPos)', nbrSamples, 1)
        name = repmat(sampleNames, 1, nbrPos)
        df = DataFrame(position=pos[:], coverage=cov[:], name=name[:])

        xTickStep = Int( 10^ceil(log10(nbrPos/2))/10 )
        xTicks = xTickStep:xTickStep:nbrPos


        covLayer = layer(df,x=:position,y=:coverage,color=:name,Geom.line)
        covPlotParams = (Scale.y_log10, Coord.cartesian(xmin=1,xmax=nbrPos,ymin=0), Guide.xticks(ticks=collect(xTicks)))

        filename = outName
        if !isempty(segmentName)
            filename = "$(filename)_$segmentName"
        end

        for f in outFormats
            push!(outFiles, _savecoverageplot(covLayer,f,joinpath(outputFolder,filename),covPlotParams...))
            push!(outSegment, segmentName)
        end
    end
    println(log, "Done.")

    outFiles, outSegment
end



# Create coverage plots for a single run
function coverageplots(syn, alignmentFolder, runName, doUpload, scriptFilename; 
                       bamCache="bam", outputFolder="plots", outFormats=[:png, :pdf], 
                       log=STDOUT, activityName = "Read Coverage Plots")
# --- Setup --------------------------------------------------------------------
    # Where to find sample .bam files. Synapse Folder ID or local folder. Synapse ID should point to "MyProject/Analysis/Alignment/Bam".
    bamPath = childpath(syn, alignmentFolder, "Bam")

    # Where to find sample .bam files. Synapse Folder ID or local folder. Synapse ID should point to "MyProject/Analysis/Alignment/Bam".
    scriptPath = childpath(syn, childpath(syn, alignmentFolder, "Scripts"), runName)

    alignLogFile = childpath(syn,scriptPath,"AlignUtils.log")
# --- Cleanup ------------------------------------------------------------------
    isdir(bamCache) || mkdir(bamCache)
    isdir(outputFolder) || mkdir(outputFolder)

# --- Make plots ---------------------------------------------------------------
    # find samples
    samplePaths, sampleNames = find_aligned(syn, bamPath, runName)
    references = referencefromlog(syn, alignLogFile, sampleNames)

    @assert !any(isempty,references)
    
    println(log, "Downloading samples")
    localPaths = map(x->localpath(syn,x,downloadLocation=bamCache,ifcollision="overwrite.local"), samplePaths)
    #sampleNamesShort = map(x->replace(x,Regex("^$(runName)_"),""), sampleNames) # remove runName_ from sample names

    filesToUpload = Vector{String}()
    fileSubFolder = Vector{String}()
    activities = Vector{Activity}() # Synapse Provenance

    # split samples by reference sequence
    for reference in unique(references)
        mask = references.==reference
        nbrSamples = countnz(mask)
        println(log, "Processing $nbrSamples samples for reference \"$reference\".")
        filesCreated, segments = _coverageplots(localPaths[mask], sampleNames[mask], outputFolder, "readcoverage_$(runName)_$reference", outFormats, log)
        for (f,s) in zip(filesCreated,segments)
            push!(filesToUpload, f)
            push!(fileSubFolder, s)
            act = Activity(activityName)
            used(act, samplePaths[mask])
            used(act, alignLogFile)
            push!(activities, act)
        end
    end

# --- Upload -------------------------------------------------------------------
    if doUpload
        basicAnalysesFolder = createchildfolder(syn, alignmentFolder, "Basic Analyses")
        runUploadFolder = createchildfolder(syn, basicAnalysesFolder, runName)
        plotUploadFolder = createchildfolder(syn, runUploadFolder, "Read Coverage")

        # upload script
        println(log, "Uploading script file.")
        scriptID = uploadiflocal(syn, scriptPath, scriptFilename, activityName=activityName)

        # upload plots
        println(log, "Uploading plots.")
        for (f, s, act) in zip(filesToUpload, fileSubFolder, activities)
            uploadFolder = isempty(s) ? plotUploadFolder : createchildfolder(syn, plotUploadFolder, s)
            executed(act, scriptID)
            uploadiflocal(syn, uploadFolder, f, act)
        end
    end

end


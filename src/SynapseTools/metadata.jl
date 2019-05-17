

function _SNumber(name)
    m = match(r"(?<=_S)\d+$", name)
    m != nothing ? m.match : ""
end


# Creates DataFrame with columns (SampleID, Run, Reference, S (if applicable), Rejected, Comment)
# parent - Bam folder
# runName - Name of run to create template for
# logFile - "AlignUtils.log", used to find reference used for each sample
function metadatatemplate(syn, runName::AbstractString; bamFolder::AbstractString="", logFile::AbstractString="", exclude::Pattern=r"Undetermined")
    isempty(bamFolder) && error("Bam folder must be specified.")

    # 1. List aligned samples
    # 2. Ignored excluded
    # 3. Create table with columns
    #   a. SampleID
    #   b. Run
    #   c. Reference (use Provenance?)
    #   d. S number (if applicable)
    #   e. Rejected
    #   f. Comment

    ids,names = listaligned(syn, bamFolder, runName::AbstractString)
    mask = .~matchany(exclude, names)
    ids,names = ids[mask],names[mask]
    map!(x->x[1:end-4], names) # remove ".bam" from end

    references = isempty(logFile) ? "" : referencefromlog(syn, logFile, names)

    # extract S numbers
    SNumbers = map(_SNumber,names)

    DataFrame(SampleID=names, Run=runName, Reference=references, S=SNumbers, Rejected="", Comment="")
end


function metadatatemplate(syn, alignmentFolder::AbstractString, runName::AbstractString; exclude::Pattern=r"Undetermined", logFileName="AlignUtils.log")
    bamFolder = getchildbyname(syn, alignmentFolder, "Bam")
    scriptsFolder = getchildbyname(syn, alignmentFolder, "Scripts")
    scriptsRunFolder = getchildbyname(syn, scriptsFolder, runName)
    logFile = getchildbyname(syn, scriptsRunFolder, logFileName)

    metadatatemplate(syn,runName,bamFolder=bamFolder,logFile=logFile,exclude=exclude)
end



# filtermetadata(metadata, filter1, filter2, ...)
# where filter is one of the following
#        :fieldname                - Check if the fieldname is a column in the metadata DataFrame. Otherwise return empty DataFrame.
#        :fieldname=>value         - Only keep samples that has the given value for the given field.
#        :fieldname=>[v1,v2...,vn] - Only keep samples that one of the given values for the given field.
#        :fieldname=>Predicate     - Only keep samples where the predicate function evaluates to true for the value of the given field name.
#        Predicate                 - Only keep samples for which Predicate returns true. Predicate acts on entire sample row.
filtermetadata(metadata::DataFrame, filter::Symbol) = filter in names(metadata) ? metadata : DataFrame()
function filtermetadata(metadata::DataFrame, filter::Pair{Symbol,T}) where {T<:AbstractArray}
    filter.first in names(metadata) || return DataFrame()
    mask = map!(x->!ismissing(x) && x in filter.second, BitArray(size(metadata,1)), metadata[filter.first])
    countnz(mask)==0 ? DataFrame() : metadata[mask,:]
end
function filtermetadata(metadata::DataFrame, filter::Pair{Symbol,T}) where {T}
    filter.first in names(metadata) || return DataFrame()
    mask = .~ismissing.(metadata[filter.first]) .& (metadata[filter.first] .== filter.second)
    countnz(mask)==0 ? DataFrame() : metadata[mask,:]
end
function filtermetadata(metadata::DataFrame, filter::Pair{Symbol,Function}) 
    filter.first in names(metadata) || return DataFrame()
    mask = map!(filter.second, BitArray(size(metadata,1)), metadata[filter.first])
    countnz(mask)==0 ? DataFrame() : metadata[mask,:]
end
function filtermetadata(metadata::DataFrame, filter::Function) 
    N = size(metadata,1)
    mask = BitArray(N)
    for i=1:N
        mask[i] = filter(metadata[i,:])
    end
    countnz(mask)==0 ? DataFrame() : metadata[mask,:]
end

function filtermetadata(metadata, args...)
    for arg in args
        metadata = filtermetadata(metadata,arg)
    end    
    metadata
end



function _loadmetadata(syn, fileID) 
    T = readtable(localpath(syn, fileID))
    T[:MetadataID] = fileID # make an extra column with a reference to the metadata file!
    T
end



getsamplemetadata(syn, file::File, args...) = filtermetadata(_loadmetadata(syn, file.id), args...)
function getsamplemetadata(syn, folder::Folder, args...) 
    #filtermetadata(_loadmetadata(syn, metadataID), args...)

    paths, names = listfiles(syn, folder.id) # TODO: update listfiles to support arguments of type Folder 
    mask = map(x->lowercase(splitext(x)[2])==".csv", names)
    paths, names = paths[mask], names[mask]

    allMetadata = Vector{DataFrame}(length(paths))
    map!(path->getsamplemetadata(syn,path,args...), allMetadata, paths)
    vcat(allMetadata...)
end

# metadataID: Synapse id, could be file or folder (should point to MyProject/Metadata).
#             If it's a folder, all .csv files in the folder will be checked.
#             filtering specified by args... (see filtermetadata())
getsamplemetadata(syn, metadataID::AbstractString, args...) = getsamplemetadata(syn, get(syn, metadataID, downloadFile=false), args...)





function appendsynapseids!(syn, metadata::DataFrame, folderID::AbstractString, fileSuffixes::AbstractArray, columnNames::AbstractArray{Symbol}; subFolder="")
    @assert length(fileSuffixes)==length(columnNames)
    N = size(metadata,1)
    synapseIDs = Vector{String}(N)

    for columnName in columnNames
        metadata[columnName] = ""
    end

    # partition samples by run, to make as few Synapse queries as possible
    runNames = metadata[:Run]
    for runName in unique(runNames)
        runInd = find(runNames.==runName) # index of samples in this run in metadata table
        sampleNames = metadata[runInd,:SampleID]

        sampleFolder = getchildbyname(syn, folderID, runName)
        isempty(subFolder) || (sampleFolder = getchildbyname(syn, sampleFolder, subFolder))

        fileIDs,filenames = listfiles(syn, sampleFolder)

        
        for (suffix,columnName) in zip(fileSuffixes,columnNames)
            fileMask = BitArray(length(filenames))
            map!( x->endswith(x,suffix), fileMask, filenames ) # find matches for suffix
            
            matchingIDs   = fileIDs[fileMask] # Synapse ids for files with matching suffix
            matchingNames = map(x->x[1:end-length(suffix)], filenames[fileMask]) # get rid of suffix from name

            matchingInd = indexin(sampleNames, matchingNames)

            matchingMask = matchingInd.!=0
            metaInd = runInd[matchingMask]
            matchingInd = matchingInd[matchingMask]

            metadata[metaInd,columnName] = matchingIDs[matchingInd]
        end

    end

    metadata
end
appendsynapseids!(syn, metadata::DataFrame, folderID::AbstractString, fileSuffix::AbstractString, columnName::Symbol; kwargs...) = appendsynapseids!(syn, metadata, folderID, [fileSuffix], [columnName]; kwargs...)


# appendswarmids!(syn, metadata::DataFrame, alignmentFolderID::AbstractString) = appendsynapseids!(syn, metadata, getchildbyname(syn,alignmentFolderID,"MutantSwarms"),".jld",:SwarmID)

# strands is one of: :both, :forward, :reverse, or an array of these
function appendswarmids!(syn, metadata::DataFrame, alignmentFolderID::AbstractString, strand::Symbol=:both)
    @assert strand in [:both, :forward, :reverse]
    if strand==:both
        appendsynapseids!(syn, metadata, getchildbyname(syn,alignmentFolderID,"MutantSwarms"), ".jld", :SwarmID)
    elseif strand==:forward
        appendsynapseids!(syn, metadata, getchildbyname(syn,alignmentFolderID,"MutantSwarms"), ".jld", :SwarmForwardID, subFolder="forward")
    else#if strand==:reverse
        appendsynapseids!(syn, metadata, getchildbyname(syn,alignmentFolderID,"MutantSwarms"), ".jld", :SwarmReverseID, subFolder="reverse")
    end
end
function appendswarmids!(syn, metadata::DataFrame, alignmentFolderID::AbstractString, strands::AbstractArray)
    for s in strands
        metadata = appendswarmids!(syn, metadata::DataFrame, alignmentFolderID::AbstractString, s)
    end
    metadata
end


function appendalignedids!(syn, metadata::DataFrame, alignmentFolderID::AbstractString)
    appendsynapseids!(syn, metadata, getchildbyname(syn,alignmentFolderID,"Bam"),
                      [".bam",".bam.bai","_consensus.fasta"],
                      [:BamID,:BaiID,:ConsensusID])
end



function downloadbymeta!(syn, metadata::DataFrame, IDColumn::Symbol, pathColumn::Symbol; cacheDir="")
    @assert haskey(metadata,IDColumn)
    
    kwargs = isempty(cacheDir) ? () :    

    ids = metadata[IDColumn]
    mask = .~ismissing.(ids)

    paths = isempty(cacheDir) ?
            localpath(syn, ids[mask]) :
            localpath(syn, ids[mask], downloadLocation=cacheDir, ifcollision="overwrite.local")

    metadata[pathColumn] = ""
    metadata[mask,pathColumn] = paths
    metadata
end

#downloadswarms!(syn, metadata::DataFrame; kwargs...) = downloadbymeta!(syn, metadata, :SwarmID, :SwarmPath; kwargs...)

# strands is one of: :both, :forward, :reverse, or an array of these
function downloadswarms!(syn, metadata::DataFrame, strand::Symbol=:both; cacheDir="", kwargs...)
    @assert strand in [:both, :forward, :reverse]
    if strand==:both
        downloadbymeta!(syn, metadata, :SwarmID, :SwarmPath; cacheDir=cacheDir, kwargs...)
    elseif strand==:forward
        isempty(cacheDir) || (cacheDir=joinpath(cacheDir,"forward"))
        downloadbymeta!(syn, metadata, :SwarmForwardID, :SwarmForwardPath; cacheDir=cacheDir, kwargs...)
    else#if strand==:reverse
        isempty(cacheDir) || (cacheDir=joinpath(cacheDir,"reverse"))
        downloadbymeta!(syn, metadata, :SwarmReverseID, :SwarmReversePath; cacheDir=cacheDir, kwargs...)
    end
end
function downloadswarms!(syn, metadata::DataFrame, strands::AbstractArray; kwargs...)
    for s in strands
        metadata = downloadswarms!(syn, metadata, s; kwargs...)
    end
    metadata
end


downloadconsensuses!(syn, metadata::DataFrame; kwargs...) = downloadbymeta!(syn, metadata, :ConsensusID, :ConsensusPath; kwargs...)
downloadbams!(syn, metadata::DataFrame; kwargs...) = downloadbymeta!(syn, metadata, :BamID, :BamPath; kwargs...)
downloadbais!(syn, metadata::DataFrame; kwargs...) = downloadbymeta!(syn, metadata, :BaiID, :BaiPath; kwargs...)
function downloadaligned!(syn, metadata::DataFrame; kwargs...) 
    downloadbams!(syn, metadata; kwargs...)
    downloadbais!(syn, metadata; kwargs...)    
    downloadconsensuses!(syn, metadata; kwargs...)
end

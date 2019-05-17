# List samples after alignment. Looks for sampleIDs by findind files called "SAMPLE_ID.bam".
# syn  - only used if path is a synapseID
# path - a local path or a synapseID
# runName - name of the run. Only used if path is a SynapseID, where it assumes runName is a subfolder. Otherwise path is used directly.
# returns paths, names (where paths are local paths or Synapse IDs)
function listaligned(syn, path::AbstractString, runName::AbstractString)
    if SynapseClient.utils.is_synapse_id(path) != nothing
        folder = getchildbyname(syn, path, runName)
        isempty(folder) && error("Could not find \"$runName\" in \"$path\".")
    else
        folder = path
    end
    paths, names = listfiles(syn, folder)

    # filter by extension
    # mask = map( x->ismatch(r".bam$",x), names )
    mask = match.(r".bam$",names) .!= nothing
    paths[mask], names[mask]
end


# get sample names and references used for all samples in logFile
function referencefromlog(syn, logFile::AbstractString)
    logFileText = readlines(localpath(syn, logFile)) # load entire text file
    map!(rstrip, logFileText, logFileText)

    # find lines of type
    # Assigned reference "RefName" to sample "SampleName".

    pattern = r"^Assigned reference \"([^\"]+)\" to sample \"([^\"]+)\".$"
    matches = map( x->match(pattern,x), logFileText )
    matches = matches[map!(x->x!=nothing, BitVector(length(matches)), matches)] # remove non-matches

    logRefs = [splitext(m.captures[1])[1] for m in matches] # reference name without file ending (.fasta)
    logSamples = [m.captures[2] for m in matches]           # sample names

    logSamples, logRefs
end


# get references used for given logFile and list of samples
function referencefromlog(syn, logFile::AbstractString, sampleNames::AbstractArray)
    logSamples, logRefs = referencefromlog(syn, logFile)

    references = Array{String}(length(sampleNames))
    for (i,name) in enumerate(sampleNames)
        ind = findfirst(logSamples, name)
        references[i] = ind!=0 ? logRefs[ind] : ""
    end
    references
end

# get references used for given logFile and a single sample
function referencefromlog(syn, logFile::AbstractString, sampleNames::AbstractString) 
    referencefromlog(syn, logFile, [sampleNames])[1]
end


function wrongreferencefromlog(syn, logFile::AbstractString)
    logFileText = readlines(localpath(syn, logFile)) # load entire text file
    map!(rstrip, logFileText, logFileText)

    # find lines of type
    # WARNING: SampleName is closer to reference "wrongRef.fasta" than to reference "RefName.fasta".

    pattern = r"^WARNING: ([^\"]+) is closer to reference \"([^\"]+)\" than to reference \"([^\"]+)\"\.$"
    matches = map( x->match(pattern,x), logFileText )
    matches = matches[map!(x->x!=nothing, BitVector(length(matches)), matches)] # remove non-matches

    logSamples = String[m.captures[1] for m in matches]     # sample names
    logWrongRefs = [splitext(m.captures[2])[1] for m in matches] # reference name without file ending (.fasta)
    logOrigRefs  = [splitext(m.captures[3])[1] for m in matches] # reference name without file ending (.fasta)

    logSamples, logWrongRefs, logOrigRefs
end

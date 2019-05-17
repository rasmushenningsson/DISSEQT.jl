


# fitnessPath: Synapse id, file or folder (should point to MyProject/Metadata).
#              If it's a folder, all .csv files in the folder will be checked.
function getsamplefitness(syn, fitnessPath)
    fitness = getsamplemetadata(syn, fitnessPath) # hijack metadata loading which is very similar
    rename!(fitness, :MetadataID, :FitnessID)
end


function appendfitness!(metadata::DataFrame, fitnessTable::DataFrame; allowInexactMatches=true)
    ids = metadata[:SampleID]
    @assert length(unique(ids)) == length(ids)
    ids = dropna(fitnessTable[:SampleID])
    @assert length(unique(ids)) == length(ids)
    @assert eltype(metadata[:Passage])<:Integer
    @assert eltype(fitnessTable[:Passage])<:Integer

    fitnessTable = copy(fitnessTable) # create copy since we will modify the table
    fitness = hcat(fitnessTable[:A],fitnessTable[:B],fitnessTable[:C])  # create Nx3 DataArray
    fitness = map(i->mean(dropna(fitness[i,:][:])), 1:size(fitness,1)); # take mean of each row after dropping NAs
    fitnessTable[:Fitness] = fitness
    fitnessTable = fitnessTable[.~isnan.(fitnessTable[:Fitness]), :] # drop measurements where 3 fitness values for a sample are NAs


    # create new columns
    metadata[:Fitness] = 0.0
    metadata[:FitnessID] = ""

    # find exact SampleID matches between metadata and fitnessTable
    matches = indexin(metadata[:SampleID], map(string,fitnessTable[:SampleID])) # converting to string converts NA->"NA" which can be handled by indexin
    hasMatches = matches.!=0
    matches = matches[hasMatches]

    metadata[hasMatches, :Fitness]   = fitnessTable[matches, :Fitness]
    metadata[hasMatches, :FitnessID] = fitnessTable[matches, :FitnessID]

    if allowInexactMatches
        fitnessTable = fitnessTable[setdiff(1:end,matches),:] # remove those fitness entries that have been used
        used = falses(size(fitnessTable,1)) # but from now on, keep track of which ones we have used rather than removing rows in the loop

        for i in find(ismissing.(metadata[:Fitness])) # for those samples we don't have a fitness value for
            if ismissing(metadata[i,:Reference]) || ismissing(metadata[i,:Mutagen]) || ismissing(metadata[i,:Dose]) || ismissing(metadata[i,:Passage])
                continue # skip samples with incomplete metadata
            end

            mask = (metadata[i,:Reference] .== fitnessTable[:Virus]) .&
                   (metadata[i,:Mutagen]   .== fitnessTable[:Mutagen]) .&
                   (metadata[i,:Dose]      .== fitnessTable[:Dose]) .&
                   (metadata[i,:Passage]   .== fitnessTable[:Passage]) .&
                   .~used

            ind = findfirst(mask) # if there are multiple matches, use the first one and leave others to match other samples
            ind==0 && continue # no match

            metadata[i,:Fitness]   = fitnessTable[ind,:Fitness]
            metadata[i,:FitnessID] = fitnessTable[ind,:FitnessID]
            used[ind] = true
        end

        # fitnessTable = fitnessTable[.~used,:]
    end

    metadata
end

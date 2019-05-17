

function id2name(syn, id)
    e = get(syn, id, downloadFile=false)
    e.name
end
function id2name(syn, id, version)
    e = get(syn, id, downloadFile=false, version=version)
    e.name
end



function retrystore(syn::Synapse, args...; kwargs...)
    delays = Int[5, 30, 60]
    i = 0 # how many times did we fail?
    while true
        try
            return store(syn, args...; kwargs...)
        catch ex
            i += 1
            if i>length(delays)
                println("Synapse store failed. Aborting.")
                rethrow(ex)
            end

            println("Synapse store failed. Waiting $(delays[i])s before trying again.")
            sleep(delays[i])
        end
    end
end




function getchildbyname(syn, parentID::AbstractString, child::AbstractString)::String
    # if 0 results, return "",
    # if 1 result, return it
    # if 2 results, error

    res = ""
    for c in getchildren(syn, parentID)
        if get(c,"name",nothing) == child
            @assert isempty(res) "Unexpected error, multiple children with the same name."
            res = c["id"]
        end
    end
    res

    # results = chunkedquery(syn, "select id from entity where entity.parentId=='$parentID' and entity.name=='$child'")
    # res = ""
    # for (i,r) in enumerate(results)
    #     i>1 && error("Unexpected error, multiple children with the same name.")
    #     res = r["entity.id"]
    # end
    # res::String
end
getchildbyname(syn, parent::AbstractEntity, child::AbstractString) = getchildbyname(syn, parent.id,child)
getchildbyname(syn, parent, args::AbstractString...) = getchildbyname(syn,getchildbyname(syn, parent, args[1]),args[2:end]...)



# Get path or SynapseID of child.
# syn        - only used if parentPath is a synapseID
# parentPath - a local path or a synapseID
# childName  - the name of a subfolder 
function childpath(syn, parentPath::AbstractString, childName::AbstractString)
	if SynapseClient.utils.is_synapse_id(parentPath) != nothing
		getchildbyname(syn, parentPath, childName)
	else
		joinpath(parentPath, childName)
	end
end



# create child folder(s) if it doesn't already exist
function createchildfolder(syn, parentID::AbstractString, childNames::AbstractString...)
    i=1
    while i≤length(childNames) # go down until we don't find a folder by that name
        childID = getchildbyname(syn, parentID, childNames[i])
        isempty(childID) && break
        parentID=childID
        i+=1
    end
    while i≤length(childNames) # create folders until end of list
        folder = Folder(childNames[i], parent=parentID)
        folder = retrystore(syn, folder)
        parentID = folder.id
        i+=1
    end
    parentID
end





# return id of uploaded file
function uploadiflocal(syn, destID::AbstractString, filePath::AbstractString, activity::Activity; fileName="")
    SynapseClient.utils.is_synapse_id(filePath) != nothing && return filePath # file is already in synapse

    isempty(fileName) && (fileName = splitdir(filePath)[2])

    file = File(path=filePath, name=fileName, parent=destID)
    file = retrystore(syn, file, activity=activity)
    file.id
end

function uploadiflocal(syn, destID::AbstractString, filePath::AbstractString; fileName="", activityName::AbstractString="", dependencies=[], exec::AbstractString="")
    activity = Activity(name=activityName)
    isempty(exec) || executed(activity, exec)
    if !(typeof(dependencies)<:Array) || !isempty(dependencies)
        used(activity, dependencies)
    end
    uploadiflocal(syn, destID, filePath, activity, fileName=fileName)
end





# Similar to readdir() but in Synapse. Returns vector of child ids as well as vector of child names.
function synapse_listfiles(syn::Synapse, parentID::AbstractString)
    ids = Vector{String}()
    names = Vector{String}()
    for c in getchildren(syn, parentID)
        push!(ids,c["id"])
        push!(names,c["name"])
    end
    ids, names

    # results = chunkedquery(syn, "select id, name from entity where entity.parentId=='$parentID'")

    # ids = Vector{String}()
    # names = Vector{String}()
    # for r in results
    #     push!(ids, r["entity.id"])
    #     push!(names, r["entity.name"])
    # end

    # ids, names
end
synapse_listfiles(syn::Synapse, parent::AbstractEntity) = synapse_listfiles(syn, parent.id)



# List files in given path. Returns both file paths and file names.
# syn  - only used if path is a synapseID
# path - a local path or a synapseID
# returns paths, names (where paths are local paths or Synapse IDs)
function listfiles(syn, path::AbstractString)
    if SynapseClient.utils.is_synapse_id(path) != nothing
        synapse_listfiles(syn, path)
    else
        names = readdir(path)
        paths = map( x->joinpath(path,x), names )
        paths, names
    end
end






# Get local path. Downloads file from Synapse if necessary.
function localpath(syn, path::AbstractString; kwargs...)
    if SynapseClient.utils.is_synapse_id(path) != nothing
        file = get(syn, path; kwargs...) # downloadLocation=cacheDir, ifcollision="overwrite.local"
        file.path
    else
        path
    end
end

localpath(syn, paths::AbstractArray; kwargs...) = map!(x->localpath(syn,x;kwargs...), Vector{String}(length(paths)), paths)






# UploadList is a vector of (uploadFolder, files, dependencies)
const UploadList = Vector{Tuple{String, Array{String}, Array{String}}}
markforupload!(list, uploadFolder::AbstractString, files::Array{T}, dependencies...) where {T<:AbstractString} =
    push!(list, (convert(String,uploadFolder), convert(Array{String},files), convert(Array{String},vcat(dependencies...))))
markforupload!(list, uploadFolder::AbstractString, file::AbstractString, dependencies...) = markforupload!(list, uploadFolder, [file], dependencies...)

function uploadfiles(syn, uploadFolder::AbstractString, uploads::UploadList, scriptFileName::AbstractString, activityName::AbstractString)
    subFolderIDs = Dict{String,String}() # "path/with/slashes" => id
    for subFolder in unique(map(x->x[1], uploads))
        subFolderIDs[subFolder] = createchildfolder(syn, uploadFolder, split(subFolder,'/',keep=false)... ) # keep=false handles subFolder="" which specifies the root (uploadFolder)
    end

    # upload script file
    scriptFile = uploadiflocal(syn, uploadFolder, scriptFileName, activityName=activityName)

    for u in uploads
        folder = subFolderIDs[u[1]]
        for f in u[2]
            uploadiflocal(syn, folder, f, activityName=activityName, dependencies=u[3], exec=scriptFile)
        end
    end
end
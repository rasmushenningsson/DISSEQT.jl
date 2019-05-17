# for convenience (use DISSEQT typealiases instead?)
const Seq = Tuple{String,String}
const Reference = Array{Seq,1}


function askforconfirmation(str::AbstractString)
    while true
        println(str, " (y/n)")
        response = strip(readline())
        lowercase(response) in ["y","yes"] && return true
        lowercase(response) in ["n","no"] && return false
    end
end


# Makes sure folder exists and is empty. 
# [Default] Asks for confirmation before removing any files.
# return false if failed
function makecleanfolder(folder::AbstractString, askForConfirmation=true)
    if isdir(folder) && !isempty(readdir(folder))
        (askForConfirmation && askforconfirmation("\"$folder\" folder is non-empty, delete all files?")) || return false
        rm(folder, recursive=true)    
    end
    isdir(folder) || mkdir(folder)
    return true
end



# expands (pattern, filename) to (pattern, filename, filepath, localfilepath)
function getreferenceinfo(syn, ref::Tuple{T,String}, referenceFolder::String) where {T}
    path = childpath(syn, referenceFolder, ref[2])
    ref[1], ref[2], path, expanduser(localpath(syn,path,downloadLocation=synapseCacheDir,ifcollision="overwrite.local"))
end
function getreferenceinfo(syn, refs::Vector, referenceFolder::String)
    [getreferenceinfo(syn,r,referenceFolder) for r in refs]
end
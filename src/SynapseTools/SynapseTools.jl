module SynapseTools


using SynapseClient
using ..DISSEQT
using Statistics
using DataFrames
using Missings
using CSV

import SynapseClient: AbstractEntity, Entity, File, Folder, Project, Activity


export
    retrystore,
    getchildbyname,
    childpath,
    createchildfolder,
    uploadiflocal,
    listfiles,
    localpath,
    UploadList,
    markforupload!,
    uploadfiles,
    listaligned,
    referencefromlog,
    wrongreferencefromlog,
    getreferencegenomes,
    metadatatemplate,
    filtermetadata,
    getsamplemetadata,
    appendsynapseids!,
    appendswarmids!,
    appendalignedids!,
    downloadbymeta!,
    downloadswarms!,
    downloadconsensuses!,
    downloadbams!,
    downloadbais!,
    downloadaligned!,
    getsamplefitness,
    appendfitness!


include("pattern.jl")
include("paths.jl")
include("samples.jl")
include("referencegenomes.jl")
include("metadata.jl")
include("fitness.jl")

end
using Distributed
nprocs()==1 && addprocs()
using SynapseClient
using DISSEQT.AlignUtils
using DISSEQT.SynapseTools




function main()
# --- Synapse login ------------------------------------------------------------
    # syn = nothing # set syn=nothing if working locally
    syn = SynapseClient.login()

# --- Setup --------------------------------------------------------------------
    # Set clean=true to rerun swarm inference. For clean=false, swarm inference will only be run if the log file is missing (i.e. no previous swarm inference was done).
    clean = false

    projectFolder = "syn11639899" # VignuzziLabPublic/Projects/FitnessLandscapes
    alignmentFolder = getchildbyname(syn, projectFolder, "Analysis", "Alignment")

    # Set uploadPath="synapseID" to upload files after running. Should point to "MyProject/Analysis/Alignment" folder. Set upload=nothing to skip uploading.
    uploadPath = alignmentFolder

    # Where to find sample .bam files. Synapse Folder ID or local folder. Synapse ID should point to "MyProject/Analysis/Alignment/Bam".
    bamPath = getchildbyname(syn, alignmentFolder, "Bam")

    # Name of the run. (Should corrsespond to a subfolder of bamPath if in Synapse.)
    runName = "H03UAAFXX"

    # Local logFile.
    logFile = "MutantSwarms.log"

    # Specify which strands to include in each run of the swarm computations. 
    strands = [:both, :forward, :reverse]

# --- Cleanup ------------------------------------------------------------------
    clean = clean || !isfile(logFile)

    if clean
        makecleanfolder("MutantSwarms") || return
        :forward in strands && mkdir("MutantSwarms/forward") # subfolder for forward strand
        :reverse in strands && mkdir("MutantSwarms/reverse") # subfolder for reverse strand
        isdir("bam") || mkdir("bam") # bam folder might not exist if we are downloading files from Synapse.
    else
        println("Skipping Swarm Inference [clean=false]")
    end


# --- Compute Mutant Swarms ----------------------------------------------------
    # find samples
    samplePaths, sampleNames = find_aligned(syn, bamPath, runName)
    
    strand2folder = Dict(:both=>"MutantSwarms",:forward=>"MutantSwarms/forward",:reverse=>"MutantSwarms/reverse")
    if clean
        # compute codon frequencies
        log = open(logFile, "w")
        for strand in strands
            outFolder = strand2folder[strand]

            println(log,"--- Computing mutant swarms for $strand strand(s) ---"); flush(log)
            computecodonfrequencies(syn, samplePaths, sampleNames, outFolder, strands=strand, log=log)
        end
        close(log)
    end


# --- Upload -------------------------------------------------------------------
    if uploadPath != nothing
        uploadswarms(syn, uploadPath, runName, @__FILE__, logFile, samplePaths, sampleNames, strands)
    end

end

main()

using Distributed
nprocs()==1 && addprocs()
using DISSEQT
using DISSEQT.AlignUtils
using SynapseClient
using DISSEQT.SynapseTools
using JLD
using Dates



function main()
# --- Synapse login ------------------------------------------------------------
    syn = SynapseClient.login()
    isdir("synapsecache") || mkdir("synapsecache")

# --- Setup --------------------------------------------------------------------
    # Set clean=true to rerun alignment. For clean=false, alignment will only be run if the log file is missing (i.e. no previous alignment was done).
    clean = false

    projectFolder = "syn11639899" # VignuzziLabPublic/Projects/FitnessLandscapes
    alignmentFolder = getchildbyname(syn, projectFolder, "Analysis", "Alignment")

    # Set uploadPath="synapseID" to upload files after running. Should point to "MyProject/Analysis/Alignment" folder. Set upload=nothing to skip uploading.
    uploadPath = alignmentFolder

    # Where to find sample .fastq files. Synapse Folder ID or local folder. Synapse ID should point to "MyProject/Raw Data/Sequencing".
    fastqPath = getchildbyname(syn, projectFolder, "Raw Data", "Sequencing")

    # Name of the run. (Should corrsespond to a subfolder of fastqPath.)
    runName = "H03UAAFXX"

    # This prefix will be added to all sample names. Normally same as runName. Set to "" if the .fastq files already have this prefix.
    namePrefix = runName

    # Where to find reference genomes. Synapse Folder ID or local folder. [Reference genomes will be uploaded if local.]
    # Synapse ID should point to "MyProject/Analysis/Alignment/ReferenceGenomes".
    referenceFolder = getchildbyname(syn, alignmentFolder, "ReferenceGenomes")

    # Rules for matching sample IDs to references.
    # The rule can either be a Regex or a function taking the sample ID and returning true/false.
    refs = [(r"-A1_", "Stop.fasta"), (r"-A2_", "Stop.fasta"), (r"-A3_", "Stop.fasta"), (r"-A4_", "Stop.fasta"), (r"-A5_", "WT.fasta"), (r"-A6_", "More.fasta"), (r"-A7_", "Less.fasta"), (r"-A8_", "Stop.fasta"), (r"-A9_", "WT.fasta"), (r"-A10_", "More.fasta"), (r"-A11_", "Less.fasta"), (r"-A12_", "Stop.fasta"), (r"-B1_", "Stop.fasta"), (r"-B2_", "Stop.fasta"), (r"-B3_", "Stop.fasta"), (r"-B4_", "Stop.fasta"), (r"-B5_", "WT.fasta"), (r"-B6_", "More.fasta"), (r"-B7_", "Less.fasta"), (r"-B8_", "Stop.fasta"), (r"-B9_", "WT.fasta"), (r"-B10_", "More.fasta"), (r"-B11_", "Less.fasta"), (r"-B12_", "Stop.fasta"), (r"-C1_", "Stop.fasta"), (r"-C2_", "Stop.fasta"), (r"-C3_", "Stop.fasta"), (r"-C4_", "Stop.fasta"), (r"-C5_", "WT.fasta"), (r"-C6_", "More.fasta"), (r"-C7_", "Less.fasta"), (r"-C8_", "Stop.fasta"), (r"-C9_", "WT.fasta"), (r"-C10_", "More.fasta"), (r"-C11_", "Less.fasta"), (r"-C12_", "Stop.fasta"), (r"-D1_", "Stop.fasta"), (r"-D2_", "Stop.fasta"), (r"-D3_", "Stop.fasta"), (r"-D4_", "Stop.fasta"), (r"-D5_", "WT.fasta"), (r"-D6_", "More.fasta"), (r"-D7_", "Less.fasta"), (r"-D8_", "Stop.fasta"), (r"-D9_", "WT.fasta"), (r"-D10_", "More.fasta"), (r"-D11_", "Less.fasta"), (r"-D12_", "Stop.fasta"), (r"-E1_", "Stop.fasta"), (r"-E2_", "Stop.fasta"), (r"-E3_", "Stop.fasta"), (r"-E4_", "Stop.fasta"), (r"-E5_", "WT.fasta"), (r"-E6_", "More.fasta"), (r"-E7_", "Less.fasta"), (r"-E8_", "Stop.fasta"), (r"-E9_", "WT.fasta"), (r"-E10_", "More.fasta"), (r"-E11_", "Less.fasta"), (r"-E12_", "Stop.fasta"), (r"-F1_", "Stop.fasta"), (r"-F2_", "Stop.fasta"), (r"-F3_", "Stop.fasta"), (r"-F4_", "Stop.fasta"), (r"-F5_", "WT.fasta"), (r"-F6_", "More.fasta"), (r"-F7_", "Less.fasta"), (r"-F8_", "Stop.fasta"), (r"-G1_", "Stop.fasta"), (r"-G2_", "Stop.fasta"), (r"-G4_", "Stop.fasta"), (r"-G5_", "WT.fasta"), (r"-G6_", "More.fasta"), (r"-G7_", "Less.fasta"), (r"-G8_", "Stop.fasta"), (r"-H1_", "Stop.fasta"), (r"-H2_", "Stop.fasta"), (r"-H3_", "Stop.fasta"), (r"-H4_", "Stop.fasta"), (r"-H5_", "WT.fasta"), (r"-H6_", "More.fasta"), (r"-H7_", "Less.fasta"), (r"-H8_", "Stop.fasta"),
            (r"Undetermined", "phiX174.fasta")] # Always keep this unless PhiX wasn't used in the sequencing.
    refs = getreferenceinfo(syn, refs, referenceFolder) # Get path/synapse id and local path info.

    # File with adapters. Synapse File ID or local file. [Uploaded if local.]
    adapters = getchildbyname(syn, alignmentFolder, "Adapters", "adapters.fa")

    # Local logFile.
    logFile = "AlignUtils.log"

    # Local file with sample info. Required for Uploading of a previously run alignment.
    sampleInfoFile = "sampleInfo.jld"

    # Minimum support for updating consensus during alignment. (Otherwise fallback to reference sequence.)
    consensusMinSupport = 100

    # Minimum support for updating consensus with indels during alignment. (Otherwise fallback to reference sequence.)
    consensusIndelMinSupport = 1000

    # Maximumum number of times a sample will be aligned with updated consensus.
    maxAlignIterations = 10

# --- Cleanup ------------------------------------------------------------------
    clean = clean || !isfile(logFile)

    if clean
        makecleanfolder("bam") || return
        makecleanfolder("temp") || return
    else
        println("Skipping Alignment [clean=false]")
    end


# --- Alignment ----------------------------------------------------------------

    samples = Sample[]
    if clean
        log = open(logFile,"w")
        isfile(sampleInfoFile) && rm(sampleInfoFile)

        start = time()
        println(log,"Starting alignment batch run at $(now())"); flush(log)
        samples = find_samples(syn, fastqPath, runName, namePrefix, log=log); flush(log)
        assign_reference!(samples,refs,log=log); flush(log)

        align_samples!(syn,samples,adapters,"bam","temp",log,consensusMinSupport=consensusMinSupport,consensusIndelMinSupport=consensusIndelMinSupport,maxAlignIterations=maxAlignIterations); flush(log)
        reference_sanity_check(samples, refs, log=log)

        duration = time()-start
        println(log,"Finished alignment batch run in $(duration)s."); flush(log)

        close(log)

        save(sampleInfoFile, "samples", samples, compress=true)
    end

# --- Upload -------------------------------------------------------------------
    if uploadPath != nothing
        if !clean # reload sample info if necessary
            try
                samples = load(sampleInfoFile, "samples")
            catch
                println("Sample info file is missing or corrupt. Please rerun alignment.")
            end
        end

        refPaths = [r[3] for r in refs]
        uploadaligned(syn, uploadPath, runName, @__FILE__, logFile, adapters, refPaths, samples)
    end


end

main()

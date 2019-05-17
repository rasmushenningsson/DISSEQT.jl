import SynapseClient: AbstractEntity, File, Folder, Activity


const synapseCacheDir = "synapsecache"



# upload all files related to an alignment run
# destID: should point to "MyProject/Analysis/Alignment" Synapse Folder
# runName: name of the run. Subfolder in scripts will be created for this run.
# alignScript: the script file used for alignment.
# adapters: adapters 
function uploadaligned(syn, destID::AbstractString, runName::AbstractString, 
					   alignScript::AbstractString, logFile::AbstractString,
					   adapters::AbstractString, refPaths::Vector, 
					   samples::Vector{Sample})
	# allFastqID = vcat(map(s->s.fastq,samples)...)
	allFastqID = unique(vcat(map(s->vcat(s.fastq,s.fastq2),samples)...))
	# all(x->SynapseClient.utils.is_synapse_id(x)!=nothing, allFastqID) || error("All sample fastq files must be in Synapse in order to setup Provenance.")


	activityName = "Aligned $runName"

	# --- Destination Alignment/Scripts/runName ---
	# create folder if it is needed
	scriptFolderID = createchildfolder(syn, destID, "Scripts")
	scriptRunFolderID = createchildfolder(syn, scriptFolderID, runName)

	# upload do_align.jl script
	alignScript = uploadiflocal(syn, scriptRunFolderID, alignScript, activityName=activityName)

	# upload adapters if needed
	adapters = uploadiflocal(syn, scriptRunFolderID, adapters, activityName=activityName,
	                         exec=alignScript)

	# upload reference genomes if needed
	for i=1:length(refPaths)
		refPaths[i] = uploadiflocal(syn, scriptRunFolderID, refPaths[i], activityName=activityName,
		                            exec=alignScript)
	end

	# upload log file
	logFile = uploadiflocal(syn, scriptRunFolderID, logFile, activityName=activityName,
	                        exec=alignScript, dependencies=vcat(adapters,refPaths,allFastqID))


	# --- Destination Alignment/Bam ---    
	bamFolderID = createchildfolder(syn, destID, "Bam")
	bamRunFolderID = createchildfolder(syn, bamFolderID, runName)
	# upload samples
	for s in samples
		dep = unique(vcat(s.fastq,s.fastq2,adapters,s.referencePath))

		# TODO: should .bai and .log be in Sample type instead of hardcoded here?
		files = [s.bam, join([s.bam,".bai"]), 
		         s.consensus, join([splitext(s.bam)[1],".log"])]

		for f in files
			uploadiflocal(syn, bamRunFolderID, f, activityName=activityName,
			              exec=alignScript, dependencies=dep)
		end
	end

end



# upload all files related to a swarm computation run
# destID: should point to "MyProject/Analysis/Alignment" Synapse Folder
# runName: name of the run. Subfolder in scripts will be created for this run.
# swarmScript: the script file used for computing swarms.
# logFile: log file.
# samplePaths: Sample Bam file paths (i.e. Synapse IDs)
# sampleNames: Sample identifiers
# strands: List of strands
# outFormat: List of output formats
function uploadswarms(syn, destID::AbstractString, runName::AbstractString, 
                      swarmScript::AbstractString, logFile::AbstractString,
                      samplePaths::Vector{String},
                      sampleNames::Vector{String},
                      strands::Vector{Symbol},
                      outFormat)
	typeof(outFormat)<:Array || (outFormat=[outFormat]) # ensure it is a vector

	# give error message if the .bam files are not already in Synapse
	# allFastqID = vcat(map(s->s.fastq,samples)...)
	all(x->SynapseClient.utils.is_synapse_id(x)!=nothing, samplePaths) || error("All sample bam files must be in Synapse in order to setup Provenance.")

	activityName = "Mutant Swarm Inference for $runName"

	# --- Destination Alignment/Scripts/runName ---
	# create folder if it is needed
	scriptFolderID = createchildfolder(syn, destID, "Scripts")
	scriptRunFolderID = createchildfolder(syn, scriptFolderID, runName)

	# upload compute_swarms.jl script
	swarmScript = uploadiflocal(syn, scriptRunFolderID, swarmScript, activityName=activityName)

	# upload log file
	logFile = uploadiflocal(syn, scriptRunFolderID, logFile, activityName=activityName,
	                        exec=swarmScript, dependencies=samplePaths)


	# --- Destination Alignment/MutantSwarms ---    
	swarmFolderID = createchildfolder(syn, destID, "MutantSwarms")
	swarmRunFolderID = createchildfolder(syn, swarmFolderID, runName)


	localStrandFolders = Vector{String}(length(strands))
	strandFolderIDs    = Vector{String}(length(strands))
	for (i,strand) in enumerate(strands)
        if strand==:both
        	localStrandFolders[i] = "MutantSwarms"
        	strandFolderIDs[i] = swarmRunFolderID
        else
        	localStrandFolders[i] = joinpath("MutantSwarms",string(strand))
        	strandFolderIDs[i] = createchildfolder(syn, swarmRunFolderID, string(strand))
        end
	end


	# upload sample swarms
	for (bam,name) in zip(samplePaths, sampleNames)
		for (localStrandFolder,strandFolderID) in zip(localStrandFolders,strandFolderIDs)
			for format in outFormat
				file = joinpath(localStrandFolder,string(name,'.',lowercase(string(format))))
				uploadiflocal(syn, strandFolderID, file, activityName=activityName,
				              exec=swarmScript, dependencies=bam)
			end
		end
	end
end

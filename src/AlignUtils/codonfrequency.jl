
# save: frequencies, positions, coverage, qualityThreshold, removeambiguous, algoritm,
#       method, newton regularization (if applicable)


function savedictJLD(fileprefix, dict)
	save("$fileprefix.jld", dict, compress=true)
end







function savedict(fileprefix, dict, outFormat=:JLD)
	if outFormat==:JLD
		savedictJLD(fileprefix, dict)
	else
		error("Unknown format: $outFormat")
	end
end


function computecodonfrequencies(samplePath::String, sampleName::String, 
                                 outFolder::String; 
                                 strands=:both, mappingQualityThreshold=30, 
                                 baseQualityThreshold=30, removeAmbiguous=true,
                                 method=:Newton, newtonRegularization=1e-6,
                                 maxIter=10000,
                                 outFormat=:JLD)
	log = IOBuffer()

	println(log, "Computing codon frequencies: ", sampleName)
	startTime = time()

	bamFile = BamFile(samplePath)
	freqs,positions,coverage = mlcodonfreqs(bamFile,strands=strands,
	                                        mappingQualityThreshold=mappingQualityThreshold,
	                                        baseQualityThreshold=baseQualityThreshold,
	                                        log=log,removeAmbiguous=removeAmbiguous,method=method,
	                                        newtonRegularization=newtonRegularization,
	                                        maxIter=maxIter)

	fileprefix = joinpath(outFolder,sampleName)
	d = Dict{String,Any}("codonFreqs"=>freqs,
	                          "positions"=>positions,
	                          "coverage"=>coverage,
	                          "segmentInfo"=>sequences(bamFile),
	                          "strands"=>string(strands),
	                          "mappingQualityThreshold" => mappingQualityThreshold,
	                          "baseQualityThreshold" => baseQualityThreshold,
	                          "removeAmbiguousBases" => removeAmbiguous,
	                          "algorithm" => "Maximum Likelihood",
	                          "solver" => string(method))
	method==:Newton && (d["newtonRegularization"] = newtonRegularization)


	typeof(outFormat) <: AbstractArray || (outFormat = [outFormat])
	for of in outFormat
		savedict(fileprefix,d,of)
	end

	duration = time()-startTime
	println(log, "Done in ", duration, "s.")
	String(take!(log))
end



function computecodonfrequencies(syn, samplePaths::Vector{String},
                                 sampleNames::Vector{String}, 
                                 outFolder::String; 
                                 log=STDOUT, bamDir="bam",
                                 kwargs...)

	# Threading modelled after the pmap implementation in http://docs.julialang.org/en/release-0.4/manual/parallel-computing/ (NB: not the same as actual pmap())
	# --------------------------------------------------------------------------
	procList = procs() # find the available processes
	n = length(samplePaths)
	i = 1
	# function to produce the next work item from the queue.
	# in this case it's just an index.
	nextidx() = (idx=i; i+=1; idx)
	@sync begin
		for p in procList
			if p != myid() || length(procList)==1
				@async begin
					while true
						idx = nextidx()
						if idx > n
							break
						end

						sPath = samplePaths[idx]
						sName = sampleNames[idx]

						println("Downloading sample $sName")
						sLocal = localpath(syn,sPath,downloadLocation=bamDir,ifcollision="overwrite.local")

						println("Computing codon frequencies for sample $sName")

						# compute codon frequencies in in worker thread
						logStr = fetch(@spawnat p computecodonfrequencies(sLocal,sName,outFolder;kwargs...))
						print(log, logStr); flush(log)
						println("Finished computing codon frequencies for sample $sName")
					end
				end
			end
		end
	end
	# --------------------------------------------------------------------------

	nothing
end


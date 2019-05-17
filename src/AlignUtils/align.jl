
const OneOrTwo{T} = Union{T,Tuple{T,T}}


function concatfastq(fastq::Array{String,1}, singleFastq::String; log=devnull)
	println(log, "--- Concatenating fastq from $(length(fastq)) files ---"); flush(log)
	err = DiskBuffer()
	err2 = DiskBuffer()

	cmdIn = pipeline(`gzip -dc $(fastq)`, stderr=openbuf(err))
	cmdOut = pipeline(`gzip`, stderr=openbuf(err2))

	ret = 0
	try
		run(pipeline(cmdIn, cmdOut, singleFastq))
	catch
		ret = 1
	end


	printiferror(log,closebuf(err))
	printiferror(log,closebuf(err2))
	ret != 0 && printiferror(log, "Fastq concatentation failed.")
	flush(log)
	ret
end



function _trimfastq(fIn, fOut, adapters::String, log)
	println(log, "--- Trimming fastq ---"); flush(log)

	err = DiskBuffer()

	ret = 0
	try
		run(pipeline(`fastq-mcf -l 16 -q 30 $fOut $adapters $fIn`, stdout=log, stderr=openbuf(err)))
	catch
		ret = 1
	end

	printiferror(log,closebuf(err))
	ret != 0 && printiferror(log, "Fastq trimming failed.")
	flush(log)
	ret
end
trimfastq(singleFastq::String, trimmedFastq::String, adapters::String; log=devnull) =
	_trimfastq(singleFastq, ("-o", trimmedFastq), adapters, log)
trimfastq(singleFastq::Tuple{String,String}, trimmedFastq::Tuple{String,String}, adapters::String; log=devnull) =
	_trimfastq(singleFastq, ("-o", trimmedFastq[1], "-o", trimmedFastq[2]), adapters, log)



# Utility function for extracting version info from tools that gives the version as part of the help text (printed to stderr)
function extractcmdversion(cmd::Cmd, readErr::Bool)
	out = DiskBuffer()

	# NB: bwa writes info to stderr
	try
		if readErr
			run(pipeline(cmd, stderr=openbuf(out)))
		else
			run(pipeline(cmd, stdout=openbuf(out)))
		end
	catch
	end
	str = closebuf(out)

	lines = split(str,'\n')
	m = match.(r"^Version:.*",lines)
	mask = m.!=nothing
	m = map(x->x.match, m[mask]) # keep matches only
	isempty(m) && return "Unknown Version"
	join(m, ", ") # well, shouldn't be more than one match, but it doesn't hurt to handle this case
end

bwaversion() = extractcmdversion(`bwa`, true)
samtoolsversion() = extractcmdversion(`samtools`, true)
fastqmcfversion() = extractcmdversion(`fastq-mcf -h`, false)



# unsortedBam and ref must be local paths
function bwaalign(trimmedFastq::OneOrTwo{String}, unsortedBam::String, ref::String; nbrThreads=4, log=devnull, keepUnmapped=true)
	println(log, "--- Aligning with bwa mem ---"); flush(log)
	out = DiskBuffer()
	err = DiskBuffer()

	# NB: bwa writes log to stderr, redirect to normal log file
	cmdAlign = pipeline(`bwa mem -t $nbrThreads $ref $trimmedFastq`, stderr=openbuf(out))

	args = keepUnmapped ? () : "-F4"
	cmdToBam = pipeline(`samtools view -bS $args -o $unsortedBam -`, stderr=openbuf(err)) # - means stdin
	
	ret = 0
	try
		run(pipeline(cmdAlign,cmdToBam,unsortedBam))
	catch
		ret = 1
	end

	printifinfo(log,closebuf(out))
	printifinfo(log,closebuf(err)) # samtools writes info to STDERR
	ret != 0 && printiferror(log, "Alignment failed.")
	flush(log)
	ret
end


function bwaindexfilescreated(reference::String)
	files = Array{String,1}()
	push!(files, reference * ".amb")
	push!(files, reference * ".ann")
	push!(files, reference * ".bwt")
	push!(files, reference * ".pac")
	push!(files, reference * ".sa")
end


# reference must be a local path
function bwaindex(reference::String; log=devnull)
	println(log, "--- Indexing with bwa index ---"); flush(log)
	out = DiskBuffer()
	out2 = DiskBuffer() # bwa index writes info to stderr...

	ret = 0
	try
		run(pipeline(`bwa index $reference`,stdout=openbuf(out),stderr=openbuf(out2)))
	catch
		ret = 1
	end

	printifinfo(log,closebuf(out))
	printifinfo(log,closebuf(out2))
	ret != 0 && printiferror(log, "BWA indexing failed.")
	flush(log)
	ret
end



function bamsort(unsortedBam::String, sortedBam::String; nbrThreads=4, log=devnull)
	println(log, "--- Sorting BAM ---"); flush(log)
	out = DiskBuffer()
	err = DiskBuffer()
	
	ret = 0
	try
		run(pipeline(`samtools sort -o $sortedBam -@$nbrThreads $unsortedBam`, stdout=openbuf(out), stderr=openbuf(err)))
	catch
		ret = 1
	end

	printifinfo(log,closebuf(out))
	printifinfo(log,closebuf(err)) # samtools writes info to STDERR
	ret != 0 && printiferror(log, "BAM sorting failed.")
	flush(log)
	ret
end

function bamindex(sortedBam::String; log=devnull)
	println(log, "--- Indexing BAM ---"); flush(log)
	err = DiskBuffer()
	
	ret = 0
	try
		run(pipeline(`samtools index $sortedBam`, stderr=openbuf(err)))
	catch
		ret = 1
	end

	printifinfo(log,closebuf(err)) # samtools writes info to STDERR
	ret != 0 && printiferror(log, "BAM indexing failed.")
	flush(log)
	ret
end


function removetempfiles(tempFiles::Array{String,1}; log=devnull)
	println(log, "--- Removing temporary files ---"); flush(log)
	for t in tempFiles
		if isfile(t)
			println(log,"Removing temporary file \"$t\"."); flush(log)
			rm(t)
		else
			printifwarning(log, "Couldn't remove temporary file \"$t\", file doesn't exist.")
		end
	end
	flush(log)
end

# setup for paired/unpaired
function singlefastq(fastqLocal::Vector{String}, prefix::String, tempFiles, log)
	length(fastqLocal)==1 && return 0, fastqLocal[1]
	singleFastq = "$(prefix)_concat_temp.fastq.gz"
	push!(tempFiles,singleFastq)
	ret = concatfastq(fastqLocal, singleFastq, log=log)
	ret, singleFastq
end
function singlefastq(fastqLocal::Vector{String}, fastqLocal2::Vector{String}, prefix::String, tempFiles, log)
	isempty(fastqLocal2) && return singlefastq(fastqLocal, prefix, tempFiles, log)
    ret1,s1 = singlefastq(fastqLocal,  "$(prefix)_1", tempFiles, log)
    ret1 != 0 && return ret1,(s1,"")
    ret2,s2 = singlefastq(fastqLocal2, "$(prefix)_2", tempFiles, log)
	ret2,(s1,s2)
end


function trimmedfastqname(singleFastq::String, prefix::String, tempFiles)
	name = "$(prefix)_trimmed_temp.fastq.gz"
	push!(tempFiles,name)
	name
end
trimmedfastqname(singleFastq::Tuple{String,String}, prefix::String, tempFiles) =
	(trimmedfastqname(singleFastq[1],"$(prefix)_1",tempFiles),
	 trimmedfastqname(singleFastq[2],"$(prefix)_2",tempFiles))


function align_sample!(sample::Sample, adapters::String, 
	                   outFolder::String, tempFolder::String; 
	                   log::IO=devnull, globalLog::IO=devnull, 
	                   maxAlignIterations::Int=5, consensusMinSupport::Int=100,
	                   consensusIndelMinSupport::Int=consensusMinSupport,
	                   keepUnmapped::Bool=true,
	                   nbrThreads::Int=4)
	tempFiles = Array{String,1}()

	# make sure fastq-files are concatenated
	ret,singleFastq = singlefastq(sample.fastqLocal,sample.fastqLocal2,joinpath(tempFolder,sample.name),tempFiles,log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Fastq concatenation failed.")
		return ret
	end
	
	trimmedFastq = trimmedfastqname(singleFastq, joinpath(tempFolder,sample.name), tempFiles)
	ret = trimfastq(singleFastq, trimmedFastq, adapters, log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Fastq trimming failed.")
		return ret	
	end


	unsortedBam = joinpath(tempFolder, "$(sample.name)_unsorted_temp.bam")
	push!(tempFiles,unsortedBam) # we will overwrite in successive iterations

	currReference = sample.referencePathLocal
	for i=1:maxAlignIterations

		# align
		ret = bwaalign(trimmedFastq, unsortedBam, currReference; nbrThreads=nbrThreads, log=log, keepUnmapped=keepUnmapped)
		if ret != 0
			removetempfiles(tempFiles, log=log)
			printiferror(globalLog, "$(sample.name): Alignment failed.")
			return ret	
		end


		# check if consensus equals reference
		println(log, "--- Computing consensus (minimum support for substitutions=$consensusMinSupport and indels=$consensusIndelMinSupport) ---")
		cons = consensus(unsortedBam, currReference, minSupport=consensusMinSupport, indelMinSupport=consensusIndelMinSupport)
		#ref  = collect(FastaReader(currReference))
		ref  = loadfasta(currReference)
		if isequal(cons, ref) 
			println(log, "Consensus equal to reference."); flush(log)
			break # done if they are equal
		end
		println(log, "Consensus not equal to reference."); flush(log)

		# otherwise save the new consensus as reference
		currReference = joinpath(tempFolder, "$(sample.name)_consensus_$i.fasta")
		push!(tempFiles,currReference)
		savefasta(currReference, cons)

		# and generate a bwa index for it...
		append!(tempFiles,bwaindexfilescreated(currReference))
		ret = bwaindex(currReference; log=log)
		if ret != 0
			removetempfiles(tempFiles, log=log)
			printiferror(globalLog, "$(sample.name): BWA indexing of consensus failed.")
			return ret	
		end



		if i==maxAlignIterations
			printiferror(log, "Consensus did not converge within $maxAlignIterations iterations.")
			printiferror(globalLog, "$(sample.name): Consensus failed to converge.")
			flush(log)
			return 1
		end
	end

	# copy consensus
	sample.consensus = joinpath(outFolder, sample.name * "_consensus.fasta")
	cp(currReference,sample.consensus,remove_destination=true)


	sample.bam = joinpath(outFolder, "$(sample.name).bam")
	ret = bamsort(unsortedBam, sample.bam; nbrThreads=nbrThreads, log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Bam sorting failed.")
		return ret	
	end

	ret = bamindex(sample.bam; log=log)
	if ret != 0
		removetempfiles(tempFiles, log=log)
		printiferror(globalLog, "$(sample.name): Bam indexing failed.")
		return ret	
	end

	removetempfiles(tempFiles, log=log)

	println(globalLog, "$(sample.name) aligned successfully.")
	return 0
end




# helper function made to be called by pmap
function align_single!(sample::Sample, adapters::String, 
                       outFolder::String, tempFolder::String,
                       maxAlignIterations::Int, consensusMinSupport::Int, 
	                   consensusIndelMinSupport::Int,
	                   keepUnmapped::Bool,
                       nbrThreads::Int)
	log = open(joinpath(outFolder, "$(sample.name).log"), "w")
	globalLog = IOBuffer()

	try 
		align_sample!(sample, adapters, outFolder, tempFolder, 
			          log=log, globalLog=globalLog, maxAlignIterations=maxAlignIterations,
			          consensusMinSupport=consensusMinSupport,
			          consensusIndelMinSupport=consensusIndelMinSupport,
		              keepUnmapped=keepUnmapped,
			          nbrThreads=nbrThreads)
	catch err
		printiferror(globalLog, string(err))
	end

	close(log)
	sample, String(take!(globalLog))
end




function align_samples!(syn, samples::Array{Sample,1}, adapters::String, 
                        outFolder::String, tempFolder::String,
                        log::IO;
                        maxAlignIterations::Int=5, consensusMinSupport::Int=1, 
                        consensusIndelMinSupport::Int=consensusMinSupport,
                        keepUnmapped::Bool=true,
                        nbrThreads::Int=4)

	println(log, "[bwa] ", bwaversion())
	println(log, "[samtools] ", samtoolsversion())
	println(log, "[fastq-mcf] ", fastqmcfversion())
	

	# Download adapters files from Synapse (if needed).
	adapters = localpath(syn, adapters, downloadLocation=synapseCacheDir, ifcollision="overwrite.local") # make sure adapters refers to a local file


	# Make sure all the references have bwa index files - before we start threading!!!
	refsPathsLocal = unique([s.referencePathLocal for s in samples])
	for r in refsPathsLocal
		if bwaindex(r) != 0
			println(log, "Indexing of reference \"$r\" failed.")
			return
		end
	end


	# Threading modelled after the pmap implementation in http://docs.julialang.org/en/release-0.4/manual/parallel-computing/ (NB: not the same as actual pmap())
	# --------------------------------------------------------------------------
	procList = procs() # find the available processes
	n = length(samples)
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

						s = samples[idx]
						println("Downloading sample $(s.name)")

						# download files if needed (main thread)
						# s.fastqLocal = map(f->localpath(syn,f,downloadLocation=synapseCacheDir,ifcollision="overwrite.local"), s.fastq)
						isempty(s.fastqLocal)  && (s.fastqLocal  = map(f->localpath(syn,f,downloadLocation=synapseCacheDir,ifcollision="overwrite.local"), s.fastq))
						isempty(s.fastqLocal2) && (s.fastqLocal2 = map(f->localpath(syn,f,downloadLocation=synapseCacheDir,ifcollision="overwrite.local"), s.fastq2))

						println("Aligning sample $(s.name)")

						# align in worker thread
						samples[idx], logStr = fetch(@spawnat p align_single!(s,adapters,outFolder,tempFolder,maxAlignIterations,consensusMinSupport,consensusIndelMinSupport,keepUnmapped,nbrThreads))
						print(log, logStr); flush(log)
						println("Finished aligning sample $(s.name)")
					end
				end
			end
		end
	end
	# --------------------------------------------------------------------------
end


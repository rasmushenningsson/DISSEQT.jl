mutable struct Sample
	name::String                # Avoid utf8-characters that could make us create weird filenames causing problems later.
	fastq::Vector{String}       # List of fastq-files belonging to this sample.
	fastq2::Vector{String}      # Mates if paired-end, otherwise empty.
	fastqLocal::Vector{String}  # Local path to fastq-files. Downloaded from Synapse if needed.
	fastqLocal2::Vector{String} # Mates if paired-end, otherwise empty.
	referenceName::String       # Filename of reference genome.
	referencePath::String       # Path to reference genome.
	referencePathLocal::String  # Local path to reference genome. Downloaded from Synapse if needed.

	bam::String			     # Output BAM file (sorted).
	consensus::String		     # Output consensus file.
end
Sample(name,fastq::Vector{String},fastq2::Vector{String},fastqLocal::Vector{String},fastqLocal2::Vector{String}) = 
	Sample(name,fastq,fastq2,fastqLocal,fastqLocal2,"","","","","")
Sample(name,fastq::Vector{String},fastq2::Vector{String}) = Sample(name,fastq,fastq2,String[],String[])
Sample(name,fastq::String,fastq2::String) = Sample(name,[fastq],[fastq2])
Sample(name,fastq::Vector{String}) = Sample(name,fastq,String[])
Sample(name,fastq::String) = Sample(name,[fastq])


# utility function for filling sample array from table (Matrix or DataFrame)
function find_samples(table, namePrefix::String="";
	                  name=nothing,
	                  fastq=nothing, fastq2=nothing,
	                  fastqLocal=nothing, fastqLocal2=nothing,
	                  separator=';',
	                  log=devnull)
	if namePrefix != "" && namePrefix[1] != '_'
		namePrefix = namePrefix * "_"
	end

	@assert name!=nothing  "name column must be specfified"
	@assert fastq!=nothing "fastq column must be specfified"

	samples = Sample[]

	for i=1:size(table,1)
		n = string(namePrefix, table[i,name])

		f1  = convert(Vector{String},split(table[i,fastq],separator,keep=false))
		f2  = fastq2 != nothing ? convert(Vector{String},split(table[i,fastq2],separator,keep=false)) : String[]
		f1L = fastqLocal  != nothing ? convert(Vector{String},split(table[i,fastqLocal], separator,keep=false)) : String[]
		f2L = fastqLocal2 != nothing ? convert(Vector{String},split(table[i,fastqLocal2],separator,keep=false)) : String[]

		push!(samples, Sample(n, f1, f2, f1L, f2L))
	end
	samples
end


# utility function that tries to identify sample names in folder were each sample might have multiple fastq files
function find_samples(syn, fastqPath::String, runName::String, namePrefix::String="", pattern::Regex=r".+(?=_L\d+_R1_\d+.fastq.gz$)"; log=devnull)
	if namePrefix != "" && namePrefix[1] != '_'
		namePrefix = namePrefix * "_"
	end

	sampleFolder = childpath(syn, fastqPath, runName)
	if isempty(sampleFolder)
		error("Could not find \"$runName\" in \"$fastqPath\".")
	end
	filePaths, fileNames = listfiles(syn, sampleFolder) # filePaths are local paths or Synapse IDs

	# find files matching pattern
	# matches = map( x->match(pattern,x), fileNames )
	# mask = falses(matches)
	# map!( x->x!=nothing, mask, matches );
	matches = match.(pattern,fileNames)
	mask = matches.!=nothing

	# log warning for non-matching files
	for f in fileNames[.~mask]
		printifwarning(log, "File \"$f\" was ignored.")
	end


	# remove files that do not match pattern
	fileNames = fileNames[mask]
	filePaths = filePaths[mask]
	matches = convert(Vector{RegexMatch}, matches[mask])

	# extract the matching part of the regex as strings
	matchingNames = Vector{String}(undef,size(matches))
	map!( x->x.match, matchingNames, matches )

	# put all files with the same match together (and add the run name to the sample name)
	samples = [ Sample("$namePrefix$u", filePaths[matchingNames.==u]) for u in unique(matchingNames) ]


	for s in samples
		n = length(s.fastq)
		println(log, "Found sample \"$(s.name)\" with $n fastq files.")
	end
	samples
end



# find paths (Synapse or local) and sample IDs (i.e. base file names) of aligned samples (.bam)
function find_aligned(syn, path::AbstractString, runName::AbstractString)
	paths, names = listaligned(syn, path, runName)
	names = map(x->x[1:end-4], names) # keep only matches and remove ".bam" ending
	paths, names
end




# To avoid extending Base.ismatch
_ismatch(r::Regex, s::AbstractString) = match(r,s)!=nothing
_ismatch(f::Function, s::AbstractString) = f(s)


# each sample name should match exactly one of the patterns ref[i][1] ∀i=1...n and will be assigned the corresponding reference name ref[i][2]
function assign_reference!(samples::Array{Sample,1}, refs::Vector; log=devnull)
	patterns = [r[1] for r in refs]
	refNames = [r[2] for r in refs]

	M = falses(length(refs))
	for (i,s) in enumerate(samples)
		broadcast!(_ismatch,M,patterns,[s.name])
		nbrMatches = sum(M)

		@assert nbrMatches>=1 "Sample \"$(s.name)\" did not match any reference"
		@assert nbrMatches<=1 "Sample \"$(s.name)\" matches references $(refNames[M])"

		index = findfirst(M)
		s.referenceName      = refNames[index]
		s.referencePath      = refs[index][3]
		s.referencePathLocal = refs[index][4]
		println(log,"Assigned reference \"$(s.referenceName)\" to sample \"$(s.name)\".")
	end
end

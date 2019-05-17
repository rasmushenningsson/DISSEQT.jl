

# reuse saving functions in codonfrequency.jl


function computenucleotidefrequencies(sample::String, 
                                      inFolder::String, outFolder::String; 
                                      strands=:both, mappingQualityThreshold=30, 
                                      baseQualityThreshold=30,
                                      method=:Newton, newtonRegularization=1e-6,
                                      outFormat=:JLD)
	log = IOBuffer()
	sampleName = splitext(splitdir(sample)[2])[1]

	println(log, "Computing nucleotide frequencies: ", sampleName)
	startTime = time()

	samplePath = joinpath(inFolder,sample)
	bamFile = BamFile(samplePath)
	freqs,positions,coverage = mlnucfreqs(bamFile,strands=strands,
	                                      mappingQualityThreshold=mappingQualityThreshold,
	                                      baseQualityThreshold=baseQualityThreshold,
	                                      log=log,method=method,
	                                      newtonRegularization=newtonRegularization)

	fileprefix = joinpath(outFolder,sampleName)
	d = Dict{String,Any}("nucleotideFreqs"=>freqs,
	                          "positions"=>positions,
	                          "coverage"=>coverage,
	                          "segmentInfo"=>sequences(bamFile),
	                          "strands"=>string(strands),
	                          "mappingQualityThreshold" => mappingQualityThreshold,
	                          "baseQualityThreshold" => baseQualityThreshold,
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


# TODO: logging
function computenucleotidefrequencies(samples::Array{String,1}, 
                                      inFolder::String, outFolder::String;
                                      log=STDOUT, kwargs...)

	sampleLogs = pmap(x->computenucleotidefrequencies(x,inFolder,outFolder;kwargs...),
	                  samples)

	for s in sampleLogs
		print(log, s)
		flush(log)
	end

	nothing
end

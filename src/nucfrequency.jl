function readcondnuc!(read::AbstractVector{Float64},nuc::Int,q::UInt8)
	p = phred2prob(q)
	read[:] = p/3*ones(4)
	read[nuc] = 1-p
	read
end



# Using the implementations in codonfrequency.jl. They work for nucleotides too.
# TODO: Refactor.




function setupproblem(observations::AbstractVector{NucQualityDict})
	@assert size(observations)==(4,) # allow trailing ones?

	# for each of the 4 nucleotides and each unique quality value in them
	# create one 4-vector with p(read|nucleotide)
	n = sum(length,observations)
	rcn = zeros(4,n)
	counts = zeros(1,n)

	nucCounts = zeros(4) # number of observations supporting each nucleotide, used for initial guess only

	i = 1
	for nuc=1:4
		d = observations[nuc]
		s = 0
		for (k,v) in d
			readcondnuc!(view(rcn,:,i),nuc,k)
			counts[1,i] = v
			nucCounts[nuc] += v
			s += v
			i += 1
		end
	end

	s = sum(nucCounts)
	theta0 = nucCounts/(s>0 ? s : 1)

	theta0,rcn,counts
end



# TODO: add logging!


function mlnucfreqs(observations::AbstractVector{NucQualityDict};
					method=:Newton, newtonRegularization=1e-6)
	theta0,rcn,counts = setupproblem(observations)
	coverage = sum(counts)

	all(theta0.==0) && return zeros(4), true, coverage # no (useful) reads - set all freqs to zero

	# early out if we can detect that there is support for a single nucleotide only?

	if method==:Newton
		theta,converged,nbrIter = newtonsolve(theta0,rcn,counts,regularization=newtonRegularization)
	elseif method==:GradientDescent
		theta,converged,nbrIter = gradientdescentsolve(theta0,rcn,counts)
	end

	#theta,obj,nbrIter
	theta, converged, coverage
end


# for a single sequence
function mlnucfreqs(m::NucQualityMap; log=STDOUT, method=:Newton, newtonRegularization=1e-6)
	len = size(m,2)
	freqs = zeros(4,len)
	coverage = zeros(Int,len)
	for i=1:len
		freqs[:,i], converged, coverage[i] = mlnucfreqs(m[:,i],method=method,newtonRegularization=newtonRegularization)
		converged || println(log, "WARNING: Solver did not converge (position=", i, ')')
	end
	freqs, 1:len, coverage
end


# solve for the whole genome
function mlnucfreqs(nqa::NucQualityAccumulator; log=STDOUT, method=:Newton, newtonRegularization=1e-6)
	N = length(nqa)
	freqs = Vector{Matrix{Float64}}(undef,N)
	positions = Vector{Vector{Int}}(undef,N)
	coverage = Vector{Vector{Int}}(undef,N)
	for (i,m) in enumerate(nqa)
		freqs[i],positions[i],coverage[i] = mlnucfreqs(m,log=log,method=method,
                                                       newtonRegularization=newtonRegularization)
	end

	freqs,positions,coverage
end


# solve for bam file
function mlnucfreqs(bamFile::BamFile; log=STDOUT, strands=:both,
                    mappingQualityThreshold=30, baseQualityThreshold=30,
                    ignoreChimeric::Bool=true,
                    method=:Newton, newtonRegularization=1e-6)
	nqa = nucqualitycount(bamFile, strands=strands, mappingQualityThreshold=mappingQualityThreshold,
	                      ignoreChimeric=ignoreChimeric)
	qualityfilter!(nqa, baseQualityThreshold)

	mlnucfreqs(nqa,log=log,method=method,newtonRegularization=newtonRegularization)
end

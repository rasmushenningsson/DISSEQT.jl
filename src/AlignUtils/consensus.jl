function sequencedist(x::Seq,y::Seq;sub::Int=1,ins::Int=10,del::Int=10)
	levenshtein(x[2],y[2],del,ins,sub) # high penalty for indels
end

function sequencedist(x::Reference,y::Reference;sub::Int=1,ins::Int=10,del::Int=10)
	@assert all(v->isequal(v[1][1],v[2][1]), zip(x, y)) "Cannot compare reference sequences with different names" 

	# sum distance
	sum(v->sequencedist(v[1],v[2],sub=sub,ins=ins,del=del), zip(x, y))
end


# elementwise comparison of strings
function hamming(x::Seq,y::Seq)
	length(x[2])!=length(y[2]) && return typemax(Int)
	s = 0
	for (cx,cy) in zip(x[2],y[2])
		s += cx!=cy
	end
	s
end

function hamming(x::Reference,y::Reference) # not used atm.
	sum(v->hamming(v[1],v[2]), zip(x,y))
end


function hasindels(x::Seq,y::Seq)
	length(x[2]) != length(y[2]) && return true # different length, must have indels
	sequencedist(x,y) != hamming(x,y) # if levenshtein distance is different from hamming distance, there are indels
end

function hasindels(x::Reference,y::Reference)
	any(v->hasindels(v[1],v[2]), zip(x,y))
end





function reference_sanity_check(sample::Sample,
                                references::Array{Reference}, 
                                referenceNames::Array{String,1}, log)
	if isempty(sample.consensus)
		printiferror(log, "$(sample.name): No consensus file specified. (Alignment failed?)")
		return
	end

	c = loadfasta(sample.consensus)
	
	d = zeros(Int,length(references))
	for (i,r) in enumerate(references)
		# allow some name mismatching since we are normally including PhiX
		if length(c) == 1 && length(r) == 1
			d[i] = sequencedist(c[1],r[1]) # allow name mismatch
		elseif length(c) != length(r)
			d[i] = typemax(Int)
		else
			d[i] = sequencedist(c,r)
		end
	end

	i = indmin(d)
	j = findfirst(referenceNames,sample.referenceName)

	if d[j] > d[i]
		ri = splitdir(referenceNames[i])[2]
		rj = splitdir(referenceNames[j])[2]
		printifwarning(log, "$(sample.name) is closer to reference \"$ri\" than to reference \"$rj\".")
	end

	# are there any indels (when compared to the best reference)?
	hasindels(c, references[i]) && printifwarning(log, "$(sample.name) has indels.")
	return
end



function reference_sanity_check(samples::Vector{Sample}, 
                                refs::Vector;
                                log=devnull)
	r = Reference[loadfasta(f[4]) for f in refs] # loaded fastas
	names = String[f[2] for f in refs]

	for i=1:length(samples)
		reference_sanity_check(samples[i], r, names, log)
	end
end




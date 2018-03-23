function readcondcodon!(read::AbstractArray{Float64,3},codon::CodonIndex, qt::CodonQualityTriplet)
	p = (phred2prob(qt[1]),phred2prob(qt[2]),phred2prob(qt[3]))

	v1 = p[3]/3*ones(4)
	v1[codon[3]] = 1-p[3]
	v2 = p[2]/3*ones(1,4)
	v2[codon[2]] = 1-p[2]
	v3 = p[1]/3*ones(1,1,4)
	v3[codon[1]] = 1-p[1]

	broadcast!(*,read,v1,v2,v3)
end




# TODO: combine compitations of objective, gradient and hessian since several steps are the same
objective(x,rcc,counts) = (-log.(x'*rcc)*counts)[1] # to scalar
grad(x,rcc,counts) = -(rcc ./ (x'*rcc)) * counts
constrainedgrad(g, normal=1/sqrt(length(g))*ones(g)) = g - dot(g,normal)*normal

function hessian!(H,x,rcc,counts)
	M = length(x)
	den = (rcc'x).^2

	for j=1:M
		for i=1:M
			H[i,j] = sum( (counts.*rcc[i,:].*rcc[j,:])./den ) # TODO: the numerator can be precomputed
		end
	end
	H
end



function kkt(x,rcc,counts)
	g = grad(x,rcc,counts)

	# we can compute what ν *should* be, by looking at any nonzero xᵢ
	i = indmax(x)
	ν = -g[i]

	λ = max.( (g+ν).*(x.==0), 0 )
	g,λ,ν
end

function kktprint(x,rcc,counts)
	g,λ,ν = kkt(x,rcc,counts)

	println( "minimum(x) = ", minimum(x))
	println( "∑ᵢxᵢ = $(sum(x))" )
	
	println("ν=$ν")

	println("λ: ", λ')
	println("minimum(λ) = ", minimum(λ), " (i=", indmin(λ), ")")
	println("λᵢxᵢ: ", (λ.*x)')
	println("(∇f + ∑ᵢλᵢ∇fᵢ + ν∇h)/||∇f||₂: ", (g-λ+ν)'/norm(g))
	g,λ,ν
end




function constrainproblem(x,rcc,active,keepMask=x.>0)
	x = x[keepMask]
	rcc = rcc[keepMask,:]
	active[active] = keepMask

	x,rcc,active
end


function gradientdescentsolve(x0,rcc,counts; maxIter::Integer=maxIter)
	rccFull = rcc

	x = x0
	active = trues(x)
	x,rcc,active = constrainproblem(x,rcc,active, x.>=0)

	obj = objective(x,rcc,counts)
	objFull = Inf # objective for the full problem, used to abort if no improvement can be found

	nbrIter = 0
	converged = false
	while true
		# sleep(0.1)

		if nbrIter<maxIter && countnz(active)>1 # move directly to KKT check if there is only one nonzero variable (no step is possible)
			nbrIter += 1

			# move active xᵢ that are zero and have a nonnegative gradient to inactive
			g = grad(x,rcc,counts)
			cg = constrainedgrad(g)

			cgmask = (cg.>=0) .& (x.<=0)
			# showall( [x.<=0 cg.>=0 cgmask] )
			if any(cgmask)
				x,rcc,active = constrainproblem(x,rcc,active,.~cgmask)
				obj = objective(x,rcc,counts)
				# println("Positive constrained gradient: constraining problem - ", countnz(active), " variables.")
				continue
			end

			ncg = norm(cg)
			if ncg > 1e-12 # check that constrained gradient is no too small
		 		cg = cg/ncg # normalize constrained gradient to handle scale of step size

				# find the maximum step length we can take (i.e. stepping further would make some xᵢ<0)
				t = x./cg # t value for hitting xᵢ=0 for all i
				t[cg.<=0] = Inf # moving away from zero
				maxStep = minimum(t)

				# println("maxStep: ", maxStep)

				# if taking the maximum size step lowers the objective
				# then solve the more constrained problem first
				xn = x - maxStep*cg
				xn[xn.<=0] = 0 # force to zero
				objn = objective(xn,rcc,counts)
				if objn < obj
					# setup more constrained problem
					xn[xn.<=1e-12] = 0 # a little bit more forcefull forcing to zero
					x,rcc,active = constrainproblem(xn,rcc,active)
					obj = objective(x,rcc,counts) # value of objective function is different depending on the number of active variables
					# println("Maxstep possible: constraining problem - ", countnz(active), " variables.")
					continue
				end
				


				# successively cut the step length in half until we find a lower value of the objective
				step = maxStep/2;
				objn = obj # to check if we have lowered it
				while true
					step<1e-12 && break # abort if step size is too small

					xn = x - step*cg
					objn = objective(xn,rcc,counts)

					if objn < obj
						# println("i: ", round(Int,log2(maxStep/step)), ", obj: ", obj, ", objn: ", objn, ", step: ", step)
						break
					end

					step = step/2
				end

				if objn < obj
					x = xn
					obj = objn
					# TODO: if any xᵢ is below threshold, constrain problem here too?
					continue
				end
			end
		end

		# could not improve constrained problem anymore
		# expand to full problem
		xC = x
		x = zeros(length(active))
		x[active] = xC
		rcc = rccFull
		active = trues(x)

		# check KKT conditions (for full problem)
		g,λ,ν = kkt(x,rcc,counts)
		# g,λ,ν = kktprint(x,rcc,counts)

		scaledGradCond = abs.(g-λ+ν) / norm(g) # move this to kkt?

		# move test to kkt!
		if all(scaledGradCond.<1e-6) && all(λ.>=-1e-6)
			# println("-- Done! ---")
			converged = true
			break # KKT conditions are satisfied - DONE!
		end

		nbrIter>=maxIter && break # Too many iterations. Didn't converge.


		# otherswise, try again with all variables now relaxed
		obj = objective(x,rcc,counts)
		if obj>=objFull-1e-12 # add some slight tolerance
			#println("WARNING: Could not improve objective further.")
			break
		end
		objFull = obj


		# TEST: relax only variables with a bad (g-λ+ν)/||g||
		newActive = (x.>0) .| (scaledGradCond.>=1e-6)
		# println("Relaxing ", countnz(newActive.!=(x.>0)), " variable(s). Objective: ", obj)
		x,rcc,active = constrainproblem(x,rcc,active, newActive)
		obj = objective(x,rcc,counts)



		# println("Relaxing all variables. Objective: ", obj)
	end

	x, converged, nbrIter
end








function newtonsolve(x0,rcc,counts; regularization=1e-6, maxIter::Integer=10000)
	rccFull = rcc

	x = x0
	active = trues(x)
	x,rcc,active = constrainproblem(x,rcc,active, x.>=0)

	obj = objective(x,rcc,counts)
	objFull = Inf # objective for the full problem, used to abort if no improvement can be found

	nbrIter = 0
	converged = false
	while true

		if nbrIter<maxIter && countnz(active)>1 # move directly to KKT check if there is only one nonzero variable (no step is possible)
			nbrIter += 1
			
	 		# setup "extended" hessian matrix [H 1; 1 0]
			g = grad(x,rcc,counts)
	 		M = length(x)
	 		H = zeros(M+1,M+1)
	 		H[1:M,M+1] = 1
	 		H[M+1,1:M] = 1
	 		hessian!(view(H,1:M,1:M),x,rcc,counts)
	 		# println("||H||=", vecnorm(H), ", cond(H)=", cond(H), ", ||H||₂=", norm(H))
	 		
	 		regularizer = norm(H)*regularization
	 		for i=1:M
	 			H[i,i] += regularizer # regularize hessian
	 			#H[i,i] += 1e5 #0.01 # regularize hessian # TODO: make this scale invariant
	 		end
	 		# println(H)


	 		Y = [-g; 0]
	 		V = H\Y

	 		dir = -V[1:M] # dir points in the direction of the positive gradient
	 		# println("μ: ", V[M+1])





			dirmask = (dir.>=0) .& (x.<=0)
			# showall( [x.<=0 cg.>=0 dirmask] )
			if any(dirmask)
				x,rcc,active = constrainproblem(x,rcc,active,.~dirmask)
				obj = objective(x,rcc,counts)
				# println("Trying to move outside boundary: constraining problem - ", countnz(active), " variables.")
				continue
			end


			# make sure the step doesn't take us beyond the border (i.e. make any xᵢ negative)
			xn = x - dir

			if any(xn.<=0)
				# println("Step size too big")

				t = x./dir
				maxT = minimum( t[dir.>0] )

				# println("||Δx||₂ = ", norm(dir))
				# println("maxT: ", maxT)


				dir = dir*maxT
				xn = x - dir
				xn[xn.<0] = 0 # force to zero

				objn = objective(xn,rcc,counts) # is the objective better at the border?
				if objn < obj
					# setup more constrained problem
					xn[xn.<=1e-12] = 0 # a little bit more forcefull forcing to zero
					x,rcc,active = constrainproblem(xn,rcc,active)
					obj = objective(x,rcc,counts) # value of objective function is different depending on the number of active variables
					# println("Maxstep possible: constraining problem - ", countnz(active), " variables.")
					continue
				end
			end




			ndir = norm(dir)
			# println("||Δx||₂ = ", ndir)
			# println("dir: ", dir')
			# println("∑Δxᵢ: ", sum(dir))


			# successively cut the step length in half until we find a lower value of the objective
			# step = maxStep/2
			step = 1

			objn = obj # make sure objn is initialized
			while true
				step*ndir<1e-12 && break # abort if step size is too small

				xn = x - step*dir
				objn = objective(xn,rcc,counts)

				if objn < obj
					# println("obj: ", obj, ", objn: ", objn, ", step: ", step*ndir)
					break
				end

				step = step/2
			end

			if objn < obj
				x = xn
				obj = objn
				# TODO: if any xᵢ is below threshold, constrain problem here too?
				continue
			end
		end

		# could not improve constrained problem anymore
		# expand to full problem
		xC = x
		x = zeros(length(active))
		x[active] = xC
		rcc = rccFull
		active = trues(x)

		# check KKT conditions (for full problem)
		g,λ,ν = kkt(x,rcc,counts)
		# g,λ,ν = kktprint(x,rcc,counts)

		scaledGradCond = abs.(g-λ+ν) / norm(g) # move this to kkt?

		# move test to kkt!
		if all(scaledGradCond.<1e-6) && all(λ.>=-1e-6)
			# println("-- Done! ---")
			converged = true
			break # KKT conditions are satisfied - DONE!
		end

		nbrIter>=maxIter && break # Too many iterations. Didn't converge.

		# otherwise, try again with all variables now relaxed
		obj = objective(x,rcc,counts)
		if obj>=objFull-1e-12 # add some slight tolerance (TODO: use relative tolerance instead)
			# println("WARNING: Could not improve objective further.")
			break
		end
		objFull = obj


		# # TEST: just relaxing one of the variables
		# nonActive = x.<=0
		# activeMask = x.>0
		# ind = 25#rand(find(nonActive))
		# activeMask[ind] = true
		# x,rcc,active = constrainproblem(x,rcc,active,activeMask)
		# obj = objective(x,rcc,counts)
		# println("Relaxing variable ", ind, ". Objective: ", obj)
		# g,λ,ν = kktprint(x,rcc,counts) # look at the constrained problem

		# TODO: try to relax variable with worst (g-λ+ν) instead of relaxing all???
		#       is "worst" the same for the full and for the relaxed problem?


		# TEST: relax only variables with a bad (g-λ+ν)/||g||
		newActive = (x.>0) .| (scaledGradCond.>=1e-6)
		# println("Relaxing ", countnz(newActive.!=(x.>0)), " variable(s). Objective: ", obj)
		x,rcc,active = constrainproblem(x,rcc,active, newActive)
		obj = objective(x,rcc,counts)


		# println("Relaxing all variables. Objective: ", obj)
	end

	x,converged,nbrIter
end




function setupproblem(observations::AbstractArray{CodonQualityDict,3})
	@assert size(observations)==(4,4,4) # allow trailing ones?

	# for each of the 64 codons and each unique quality triplet in them
	# create one 4x4x4 matrix with p(read|codon)
	n = sum(length,observations)
	rcc = zeros(4,4,4,n)
	counts = zeros(n)


	i = 1
	for c1=1:4
		for c2=1:4
			for c3=1:4
				d = observations[c3,c2,c1]
				# s = 0
				for (k,v) in d
					#readcondcodon!(rcc[:,:,:,i],(c1,c2,c3),k)
					readcondcodon!(view(rcc,:,:,:,i),(c1,c2,c3),k)
					counts[i] = v
					# s += v
					i += 1
				end
			end
		end
	end

	# TODO: avoid this reshape by forming it in the right way from the start
	rcc = reshape(rcc,64,n)
	# counts = reshape(counts,1,n) # keep as vector

	nbrComplete = map(nbrcompletecodons, observations)

	# initial guess
	s = sum(nbrComplete)
	theta0 = nbrComplete[:]/(s>0?s:1)

	theta0,rcc,counts
end



# TODO: add logging!


# solve for a single codon
function mlcodonfreqs(observations::AbstractArray{CodonQualityDict,3}; 
					  method=:Newton, newtonRegularization=1e-6, maxIter::Integer=10000)
	theta0,rcc,counts = setupproblem(observations)
	coverage = sum(counts)

	all(theta0.==0) && return zeros(4,4,4), true, coverage, 0 # no (useful) reads - set all freqs to zero

	# early out if we can detect that there is support for a single codon only?

	if method==:Newton
		theta,converged,nbrIter = newtonsolve(theta0,rcc,counts,regularization=newtonRegularization,maxIter=maxIter)
	elseif method==:GradientDescent
		theta,converged,nbrIter = gradientdescentsolve(theta0,rcc,counts,maxIter=maxIter)
	end

	reshape(theta,4,4,4), converged, coverage, nbrIter
end


# for a single sequence
function mlcodonfreqs(m::CodonQualityMap; log=STDOUT, method=:Newton, newtonRegularization=1e-6, maxIter::Integer=10000)
	len = size(m,4)
	freqs = zeros(4,4,4,len)
	coverage = zeros(Int,len)
	for i=1:len
		freqs[:,:,:,i], converged, coverage[i], nbrIter = mlcodonfreqs(m[:,:,:,i],method=method,newtonRegularization=newtonRegularization,maxIter=maxIter)
		converged || println(log, "WARNING: Solver did not converge (position=", i, ", nbrIter=", nbrIter, ')')
	end
	freqs, 1:len, coverage
end


# solve for the whole genome
function mlcodonfreqs(cqa::CodonQualityAccumulator; log=STDOUT, method=:Newton, newtonRegularization=1e-6, maxIter::Integer=10000)
	# Array{Float64,4}[mlcodonfreqs(m,log=log,method=method,
	# 				 newtonRegularization=newtonRegularization) for m in cqa]
	N = length(cqa)
	freqs = Array{Array{Float64,4},1}(N)
	positions = Array{Array{Int,1},1}(N)
	coverage = Array{Array{Int,1},1}(N)
	for (i,m) in enumerate(cqa)
		freqs[i],positions[i],coverage[i] = mlcodonfreqs(m,log=log,method=method,
	 				                                     newtonRegularization=newtonRegularization,
	 				                                     maxIter=maxIter)
	end

	freqs,positions,coverage
end


# solve for bam file
function mlcodonfreqs(bamFile::BamFile; log=STDOUT, strands=:both, 
					  mappingQualityThreshold=30, baseQualityThreshold=30,
					  ignoreChimeric::Bool=true,
					  removeAmbiguous=true, method=:Newton, newtonRegularization=1e-6,
					  maxIter::Integer=10000)
	cqa = codonqualitycount(bamFile, strands=strands, mappingQualityThreshold=mappingQualityThreshold,
							ignoreChimeric=ignoreChimeric)
	qualityfilter!(cqa, baseQualityThreshold)
	removeAmbiguous && removeambiguous!(cqa)

	mlcodonfreqs(cqa,log=log,method=method,newtonRegularization=newtonRegularization,maxIter=maxIter)
end

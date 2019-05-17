struct BamRegion
	bamFile::BamFile
	refID::Int # 0-based index, since that is what each read has stored
	region::UnitRange{Int}
end
#BamRegion{T<:Integer}(f::BamFile, refID::Int, region::UnitRange{T})

BamRegion(f::BamFile, refID::Int) = BamRegion(f,refID,1:sequences(f)[refID+1][2])
BamRegion(f::BamFile, refName::AbstractString) = BamRegion(f,refID(f,refName))
BamRegion(f::BamFile, refName::AbstractString, region::UnitRange) = BamRegion(f,refID(f,refName), region)

region(f::BamFile, args...) = BamRegion(f, args...)


function ChunkQueue(reg::BamRegion)
	@assert hasindex(reg.bamFile) "Missing bai index for $(filename(reg.bamFile))."
	ChunkQueue(index(reg.bamFile), reg.refID, reg.region)
end



# TODO: type or immutable?
mutable struct BamRegionState
	refID::Int
	region::UnitRange{Int}
	cq::ChunkQueue
	currChunk::Chunk
	beginOffs::VirtualOffset # linear index offset of the 16kbp window containing start of region
	s::BGZFBlockStream
	currRead::BamRead
	nextRead::BamRead
	isdone::Bool # useful since there are many indicators that we are done
end


function BamRegionState(reg::BamRegion)
	linearIndex = index(reg.bamFile).refs[reg.refID+1].linearIndex
	window = ((reg.region[1]-1)>>14)+1
	#beginOffs = linearIndex[window]
	beginOffs = window<=length(linearIndex) ? linearIndex[window] : VirtualOffset(0xffffffffffffffff)

	s = BamRegionState(reg.refID,
	                   reg.region,
	                   ChunkQueue(reg),
	                   Chunk(VirtualOffset(0),VirtualOffset(0)),
	                   beginOffs,
	                   BGZFBlockStream(filename(reg.bamFile)),
	                   BamRead(reg.bamFile),
	                   BamRead(reg.bamFile),
	                   false)
	# put in valid state such that we can call next() and done()

	nextimpl(reg.bamFile, s)
	s

	# s.done = isempty(s.cq)
	# s.done && return s # no chunks - done
	# nextimpl(s)
	# s
end




# start(reg::BamRegion) = BamRegionState(reg)

function nextimpl(bamFile::BamFile, s::BamRegionState)
	
	# println("--- nextimpl ---")

	while true

		# println("ftell: ", position(s.s))
		# println("currChunk end: ", s.currChunk.e)

		# move to a new chunk if needed
		if position(s.s) >= s.currChunk.e # after current chunk?
			while true
				isempty(s.cq) && (s.isdone=true; return)
				s.currChunk = shift!(s.cq) # get next chunk
				position(s.s) >= s.currChunk.e && continue # skip if it doesn't end after last chunk?

				# println("currChunk: ", s.currChunk)
				# println("linear index offset: ", s.beginOffs)
				
				#if true # isinteresting(currChunk)
				if s.currChunk.e > s.beginOffs # use linear index to check if chunk could be of interest
					# println("Using chunk: ", s.currChunk)

					# seek if next chunk starts after current
					position(s.s) < s.currChunk.s && seek(s.s, s.currChunk.s)
					break
				else
					# println("Skipping chunk: ", s.currChunk)
				end
			end
		end

		# println("block size: ", length(s.s.block))
		eof(s.s) && (isdone=true; return)


		# we are now in a chunk with interesting reads
		# println("nextread!")
		#nextread!(s.s,s.nextRead)
		trynextread!(s.s,s.nextRead) || (s.isdone=true; return)
		# println("after nextread!")

		# read starts beyond region?
		# println("pnext: ", pos(s.nextRead), ", region end: ", s.region[end])
		pos(s.nextRead) > s.region[end] && (s.isdone=true; return)
			
		# # read ends before region?
		# println("pendnext: ", endpos(s.nextRead), ", region start: ", s.region[1])
		# endpos(s.nextRead) < s.region[1] && continue

		# read ends before region?
		# endpos is much more expensive to compute, only do it for reads where it is needed!
		pos(s.nextRead) < s.region[1] && endpos(s.nextRead) < s.region[1] && continue


		break
	end
end


# function next(reg::BamRegion, s::BamRegionState)
# 	s.currRead,s.nextRead = s.nextRead,s.currRead # swap
# 	nextimpl(reg.bamFile, s)
# 	(s.currRead, s)
# end

# done(::BamRegion, s::BamRegionState) = s.isdone


function iterate(reg::BamRegion, s::BamRegionState=BamRegionState(reg))
	s.isdone && return nothing

	s.currRead,s.nextRead = s.nextRead,s.currRead # swap
	nextimpl(reg.bamFile, s)
	(s.currRead, s)
end


eltype(::Type{BamRegion}) = BamRead



# immutable BamState
# 	io::IO
# 	zio::IO
# 	read::BamRead
# end

# BamState(f::BamFile, io::IO, zio::IO) = BamState(io,zio,BamRead(f))
# function BamState(f::BamFile, io::IO) 
# 	zio = ZlibInflateInputStream(io,reset_on_end=true)
# 	skip_bam_header(zio) # just skip header for now
# 	BamState(f, io, zio)
# end
# BamState(f::BamFile) = BamState(f,open(f.filename))



# Base.start(f::BamFile) = BamState(f)
# Base.next(::BamFile, s::BamState) = (nextread!(s.zio, s.read), s)
# Base.done(::BamFile, s::BamState) = eof(s.zio)
# Base.eltype(::Type{BamFile}) = BamRead

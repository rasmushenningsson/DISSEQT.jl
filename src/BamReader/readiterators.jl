struct BamReadChunk
	e::CIGARElement
	readPos::Int
	refPos::Int
end

BamReadChunk() = BamReadChunk(0,0,0)


op(c::BamReadChunk) = op(c.e)
op_char(c::BamReadChunk) = op_char(c.e)
len(c::BamReadChunk) = len(c.e)

readpos(c::BamReadChunk) = c.readPos
refpos(c::BamReadChunk) = c.refPos


print(io::Main.Base.IO, c::BamReadChunk) = print(io, "Chunk: ", c.e, ", readPos: ", c.readPos, ", refPos: ", c.refPos)




struct BamReadChunkState
	cig::CIGAR
	cigState::Int
	readPos::Int
	refPos::Int
end
BamReadChunkState(cig::CIGAR) = BamReadChunkState(cig,1,1,pos(r))

# function Base.start(r::BamRead)
# 	cig = cigar(r)
# 	BamReadChunkState(cig,start(cig),1,pos(r))
# end
# function Base.next(::BamRead, s::BamReadChunkState)
# 	(e,cigIter) = next(s.cig,s.cigIter)
# 	readPos = s.readPos
# 	refPos = s.refPos

# 	o = op(e)
# 	l = len(e)
# 	if o==CIGAR_M || o==CIGAR_S || o==CIGAR_I || o==CIGAR_Eq || o==CIGAR_X
# 		readPos += l
# 	end

# 	if o==CIGAR_M || o == CIGAR_D || o==CIGAR_N || o==CIGAR_Eq || o==CIGAR_X
# 		refPos += l
# 	end
# 	# Do nothing for CIGAR_H and CIGAR_P

# 	# TODO: check for unsupported/unknown CIGAR ops

# 	(BamReadChunk(e,s.readPos,s.refPos),BamReadChunkState(s.cig,cigIter,readPos,refPos))
# end
# Base.done(::BamRead, s::BamReadChunkState) = done(s.cig,s.cigIter)


function iterate(r::BamRead, s::BamReadChunkState=BamReadChunkState(cigar(r)))
	cigIter = iterate(s.cig, s.cigState)
	cigIter == nothing && return nothing
	(e,cigState) = cigIter

	readPos = s.readPos
	refPos = s.refPos

	o = op(e)
	l = len(e)
	if o==CIGAR_M || o==CIGAR_S || o==CIGAR_I || o==CIGAR_Eq || o==CIGAR_X
		readPos += l
	end

	if o==CIGAR_M || o == CIGAR_D || o==CIGAR_N || o==CIGAR_Eq || o==CIGAR_X
		refPos += l
	end
	# Do nothing for CIGAR_H and CIGAR_P

	# TODO: check for unsupported/unknown CIGAR ops

	(BamReadChunk(e,s.readPos,s.refPos),BamReadChunkState(s.cig,cigState,readPos,refPos))
end

eltype(::Type{BamRead}) = BamReadChunk



function readrange(c::BamReadChunk)
	p = readpos(c)
	l = len(c)
	p:p+l-1
end

seq(r::BamRead,c::BamReadChunk) = view(seq(r),readrange(c))
qual(r::BamRead,c::BamReadChunk) = view(qual(r),readrange(c))



mutable struct BGZFBlockStream
	io::IO # underlying compressed data
	block::Array{UInt8,1} # current BGZF block
	#zio # TODO: set type?
	zio::BufferedInputStream{Libz.Source{:inflate,BufferedInputStream{BufferedStreams.EmptyStream}}}

	# TODO: keep track of current virtual offset?
	pos::VirtualOffset
end

function BGZFBlockStream(io::IO, block::Array{UInt8,1}=Array{UInt8,1}())
	block = readblock(io,block)
	BGZFBlockStream(io,block,ZlibInflateInputStream(block),VirtualOffset(0))
end
BGZFBlockStream(filename::AbstractString) = BGZFBlockStream(open(filename))



function readblock(io::IO, block::Array{UInt8,1})
	pos = position(io)
	header = Array{UInt8,1}(12+6)
	readbytes!(io,header,12+6)

	# println("magic: ", header[1:4])
	@assert isequal(header[1:4], UInt8[31,139,8,4])

	bSize = (UInt16(header[12+4+1])) | (UInt16(header[12+4+2])<<8)	

	seek(io,pos) # TODO: avoid this seek call
	block = resize!(block,bSize+1)
	readbytes!(io,block,bSize+1)

	# println("bSize+1=", bSize+1)

	block
end


position(s::BGZFBlockStream) = s.pos
eof(s::BGZFBlockStream) = length(s.block) == 28 # assumes the only empty block is at the end

function seek(s::BGZFBlockStream, offs::VirtualOffset)
	# TODO: check if seek is necessary
	seek(s.io, coffset(offs))

	s.block = readblock(s.io, s.block)
	s.zio = ZlibInflateInputStream(s.block)

	# seek in zio (how do I do this in a good way?)
	u = uoffset(offs)
	u > 0 && readbytes(s.zio,u)
	s.pos = offs
	return
end
eofblock(s::BGZFBlockStream) = eof(s.zio)


# assumes eof(s.zio) is false
# TODO: allow reads spanning multiple blocks
function trynextread!(s::BGZFBlockStream,r::BamRead)
	# println("s.pos=", s.pos)
	# println("eofblock(s)=", eofblock(s))
	eofblock(s) && seek(s,s.pos) # initialize next block if needed

	eof(s) && return false

	# nextread!(s.zio,r)
	# s.pos += block_size(r) + 4 # total size of read is block_size + the int32 containing the block_size value
	# # if we reached the end of the block, set virtual offset to start of next block
	# eofblock(s) && (s.pos = VirtualOffset(coffset(s.pos)+length(s.block),0))

	nextread!(s,r)

	# println("after read: s.pos=", s.pos)

	
	true
end


# TODO: make these support ints that are stored in the intersection between two blocks?
function read(s::BGZFBlockStream,::Type{Int32})
	eofblock(s) && seek(s,s.pos) # initialize next block if needed
	v = read(s.zio,Int32)
	s.pos += 4
	eofblock(s) && (s.pos = VirtualOffset(coffset(s.pos)+length(s.block),0))

	# println(s.pos)

	v
end
function read(s::BGZFBlockStream,::Type{UInt32})
	eofblock(s) && seek(s,s.pos) # initialize next block if needed
	v = read(s.zio,UInt32)
	s.pos += 4
	eofblock(s) && (s.pos = VirtualOffset(coffset(s.pos)+length(s.block),0))

	# println(s.pos)

	v
end


function readbytes!(s::BGZFBlockStream,b::Vector{UInt8},nb)
	eofblock(s) && seek(s,s.pos) # initialize next block if needed


	# println("nb: ",nb)


 	nread = readbytes!(s.zio,b,nb)
 	nb = nb-nread
	s.pos += nread
	eofblock(s) && (s.pos = VirtualOffset(coffset(s.pos)+length(s.block),0))
	# println(s.pos)

	# println("nread: ",nread)
	# println("s.pos: ",s.pos)



 	byteoffs = nread+1
 	totalRead = nread
 	while nb>0
 		eofblock(s) && seek(s,s.pos) # initialize next block if needed

 		# temporary buffer since readbytes! doesn't have API for reading into a buffer with an offset
	 	tempbuf = readbytes(s.zio,nb)
	 	nread = length(tempbuf)

		# println("nb: ",nb)
		# println("nread: ",nread)


	 	b[ byteoffs:byteoffs+nread-1 ] = tempbuf

	 	nb = nb-nread
		s.pos += nread
		eofblock(s) && (s.pos = VirtualOffset(coffset(s.pos)+length(s.block),0))
		# println(s.pos)

	 	byteoffs += nread
	 	totalRead += nread
 	end

 	totalRead
end




# # overwrite read with data from stream
# function nextread!(s::BGZFBlockStream,r::BamRead)
# 	r.block_size = read(io, Int32)
# 	r.refID      = read(io, Int32)
# 	r.pos        = read(io, Int32)
# 	r.bin_mq_nl  = read(io, UInt32)
# 	r.flag_nc    = read(io, UInt32)
# 	r.l_seq      = read(io, Int32)
# 	r.next_refID = read(io, Int32)
# 	r.next_pos   = read(io, Int32)
# 	r.tlen       = read(io, Int32)

# 	reqsize = Int(r.block_size)-8*4
# 	#reqsize > length(r.buf) && error("Read buffer too small!") # we could expand it, but if the file is corrupt, we could allocate huge things here...
# 	@assert reqsize>=0 "Corrupt Bam Read, block_size too small"

# 	# expand read buffer within reasonable limits
# 	# println("pointer: ", pointer(r.buf), ", size: ", length(r.buf))
# 	if reqsize>length(r.buf)
# 		sz = nextpow2(reqsize)
# 		@assert sz<=MAX_READ_SIZE "Read does not fit in maximum read buffer size ($MAX_READ_SIZE bytes)."
# 		r.buf = Array{UInt8,1}(sz)
# 	end

# 	println(reqsize)

# 	n = readbytes!(io,r.buf,reqsize)
# 	n!=reqsize && error("Unexpected end of stream when loading read.")
# 	compute_offsets!(r)
# 	r
# end

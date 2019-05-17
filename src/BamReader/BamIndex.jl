
# Input, 1-based, closed interval
# returned bin index is 0-based!!!
function reg2bin(s::Int, e::Int)
	@assert s>0 && e>0 "Region coordinates must be positive integers."
	@assert s<=2^29 && e<=2^29 "Region coordinates must be <= 2^29."

	s = s-1 # 1-based to 0-based
	e = e-1 # 1-based to 0-based
	# Original implementation is for closed intervals (if we skip end--), so e is correct.

	s>>14 == e>>14 && return div((1<<15)-1,7) + (s>>14)
	s>>17 == e>>17 && return div((1<<12)-1,7) + (s>>17)
	s>>20 == e>>20 && return div( (1<<9)-1,7) + (s>>20)
	s>>23 == e>>23 && return div( (1<<6)-1,7) + (s>>23)
	s>>26 == e>>26 && return div( (1<<3)-1,7) + (s>>26)
	0

	# from BAM specification pdf
	# --end;
	# if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	# if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	# if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
	# if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
	# if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
	# return 0;
end


# Input, 1-based, closed interval
# returned bin indices are 0-based!!!
function reg2bins(s::Int, e::Int)
	@assert s>0 && e>0 "Region coordinates must be positive integers."
	@assert s<=2^29 && e<=2^29 "Region coordinates must be <= 2^29."


	s = s-1 # 1-based to 0-based
	e = e-1 # 1-based to 0-based
	# Original implementation is for closed intervals (if we skip end--), so e is correct.

	bins = Array{Int,1}()
	push!(bins, 0)
	append!(bins,    1 + (s>>26) :     1 + (e>>26))
	append!(bins,    9 + (s>>23) :     9 + (e>>23))
	append!(bins,   73 + (s>>20) :    73 + (e>>20))
	append!(bins,  585 + (s>>17) :   585 + (e>>17))
	append!(bins, 4681 + (s>>14) :  4681 + (e>>14))

	# from BAM specification pdf
	#const MAX_BIN = div((1<<18)-1,7)
	# int i = 0, k;
	# --end;
	# list[i++] = 0;
	# for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list[i++] = k;
	# for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list[i++] = k;
	# for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list[i++] = k;
	# for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list[i++] = k;
	# for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	# return i;
end

# Input, 1-based, closed interval
# returned bin indices are 0-based!!!
function reg2binranges(s::Int, e::Int)
	@assert s>0 && e>0 "Region coordinates must be positive integers."
	@assert s<=2^29 && e<=2^29 "Region coordinates must be <= 2^29."


	s = s-1 # 1-based to 0-based
	e = e-1 # 1-based to 0-based
	# Original implementation is for closed intervals (if we skip end--), so e is correct.

	UnitRange{Int}[             0 :               0,
	                  1 + (s>>26) :     1 + (e>>26),
	                  9 + (s>>23) :     9 + (e>>23),
	                 73 + (s>>20) :    73 + (e>>20),
	                585 + (s>>17) :   585 + (e>>17),
	               4681 + (s>>14) :  4681 + (e>>14)]
end



struct VirtualOffset
	o::UInt64
end
function VirtualOffset(coffset::Integer, uoffset::Integer)
	@assert 0<=coffset<(UInt64(2)^48-1) && 0<=uoffset<2^16
	VirtualOffset(UInt64(coffset)<<16 | uoffset)
end

coffset(v::VirtualOffset) = v.o>>16
uoffset(v::VirtualOffset) = v.o & (1<<16-1)

# allow adding to virual offset as long as we stay within the same block
function +(v::VirtualOffset, i::Integer)
	@assert i>=0
	@assert uoffset(v)+i < (1<<16)
	VirtualOffset( v.o + i )
end

isless(u::VirtualOffset, v::VirtualOffset) = u.o<v.o


function print(io::IO, v::VirtualOffset)
	# print(io, "0x", hex(v.o,16))
	print(io, "c=0x", hex(coffset(v),12), ":u=0x", hex(uoffset(v),4))
end
show(io::IO, v::VirtualOffset) = print(io,v)



struct Chunk
	s::VirtualOffset # start
	e::VirtualOffset # end
end
function Chunk(is::IO)
	s = VirtualOffset(read(is,UInt64))
	e = VirtualOffset(read(is,UInt64))
	Chunk(s,e)
end

# define a total order of chunks
isless(a::Chunk, b::Chunk) = a.s<b.s || (a.s==b.s && a.e<b.e)

function print(io::IO, chunk::Chunk)
	#print(io, "Chunk(0x", hex(chunk.s,16), "-0x", hex(chunk.e,16), ")")
	print(io, "Chunk(", chunk.s, "-", chunk.e, ")")
end
show(io::IO, chunk::Chunk) = print(io,chunk)


mutable struct Bin
	bin::UInt # distinct bin - TODO: remove? store as key in Bin Dict only.
	chunks::Vector{Chunk}
end
function Bin(is::IO)
	bin = UInt( read(is,UInt32) )

	nChunks = Int(read(is,Int32))
	chunks = Vector{Chunk}(undef,nChunks)

	for i=1:nChunks
		chunks[i] = Chunk(is)
	end
	
	Bin(bin,chunks)
end

function print(io::IO, bin::Bin)
	print(io, "Bin(", Int64(bin.bin), ", ", bin.chunks, ")")
end
show(io::IO, bin::Bin) = print(io,bin)


# BamIndex for a single reference
mutable struct BamIndexRef
	#bins::Array{Bin,1} # TODO: change to Dict mapping distinct bin to chunk list?
	bins::Dict{Int,Bin}
	linearIndex::Array{VirtualOffset,1}
end

function BamIndexRef(is::IO)
	nBins = Int(read(is,Int32))
	#bins = Array{Bin,1}(nBins)
	bins = Dict{Int,Bin}()

	for i=1:nBins
		#bins[i] = Bin(is)
		b = Bin(is)
		bins[b.bin] = b
	end

	nIntervals = Int(read(is,Int32))
	linearIndex = reinterpret(VirtualOffset, read!(is,Vector{UInt64}(undef,nIntervals)))

	BamIndexRef(bins,linearIndex)
end


mutable struct BamIndex
	refs::Array{BamIndexRef,1}
	nbrNoCooord::UInt64
end
BamIndex() = BamIndex(Array{BamIndexRef,1}(),0) # empty bam index

isempty(index::BamIndex) = length(index.refs)==0 && index.nbrNoCooord==0

function BamIndex(is::IO)
	magic = UInt32('B') | UInt32('A')<<8 | UInt32('I')<<16 | UInt32(1)<<24
	magicRead = read(is,UInt32)

	@assert magicRead == magic "Corrupt BAI file, magic string doesn't match."

	nRef = Int(read(is,Int32))
	refs = Vector{BamIndexRef}(undef,nRef)

	for i=1:nRef
		refs[i] = BamIndexRef(is)
	end

	nbrNoCooord = UInt64(0)
	if !eof(is)
		# handle n_no_coor if it is available
		nbrNoCooord = read(is,UInt64)
	end


	BamIndex(refs,0)
end
function BamIndex(filename::AbstractString)
	is = open(filename)
	bi = BamIndex(is)
	close(is)
	bi
end


function tryloadindex(bamFilename::AbstractString)
	bai = bamFilename * ".bai" # try to load "filename.bam.bai"
	!isfile(bai) && (bai = splitext(bamFilename)[1] * ".bai") # try to load "filename.bai"
	!isfile(bai) && return BamIndex() # neither found - return empty index
	BamIndex(bai) # load index
end







mutable struct ChunkQueueSingle
	index::BamIndexRef
	binRange::UnitRange{Int}
	currBin::Int
	currChunk::Int
end

function next!(q::ChunkQueueSingle)
	q.currChunk += 1

	#while !isempty(q) && q.currChunk>length(q.index.bins[q.currBin].chunks)
	while !isempty(q) && 
	      (!haskey(q.index.bins,q.binRange[q.currBin]) || 
	       q.currChunk>length(q.index.bins[q.binRange[q.currBin]].chunks))
		q.currBin += 1
		q.currChunk = 1
	end
	q
end


function ChunkQueueSingle(index::BamIndexRef, binRange::UnitRange{Int})
	q = ChunkQueueSingle(index,binRange,1,0) # start before first
	next!(q) # move to next
end


isempty(q::ChunkQueueSingle) = q.currBin>length(q.binRange)

function shift!(q::ChunkQueueSingle)
	ret = q.index.bins[q.binRange[q.currBin]].chunks[q.currChunk] # TODO: avoid index[currBin] every iteration, unnecessary lookup.
	next!(q)
	ret
end



mutable struct ChunkQueue
	queues::Array{ChunkQueueSingle,1}
	pq::PriorityQueue{Int,Chunk} # bin level -> Chunk (the bin level is the index in array of single queues)

	function ChunkQueue(queues::Array{ChunkQueueSingle,1}, pq::PriorityQueue{Int,Chunk})
		@assert length(queues)==6 "Incorrect number of bin levels."
		for i=1:length(queues)
			isempty(queues[i]) || (pq[i] = shift!(queues[i]))
		end
		new(queues,pq)
	end
end

ChunkQueue(queues::Array{ChunkQueueSingle,1}) = ChunkQueue(queues,PriorityQueue(Int,Chunk))
ChunkQueue(index::BamIndexRef, binRanges::Array{UnitRange{Int},1}) = ChunkQueue(map(x->ChunkQueueSingle(index,x), binRanges))
ChunkQueue(index::BamIndexRef, region::UnitRange{Int}) = ChunkQueue(index,reg2binranges(region[1], region[end]))
ChunkQueue(index::BamIndex, refID::Int, region::UnitRange{Int}) = ChunkQueue(index.refs[refID+1], reg2binranges(region[1], region[end]))


isempty(q::ChunkQueue) = isempty(q.pq)

function shift!(q::ChunkQueue)
	(binLevel, chunk) = peek(q.pq)
	dequeue!(q.pq)

	# println("binLevel: ", binLevel, ", chunk: ", chunk)

	# add next chunk from this queue (if there is one)
	isempty(q.queues[binLevel]) || (q.pq[binLevel] = shift!(q.queues[binLevel]))

	chunk
end
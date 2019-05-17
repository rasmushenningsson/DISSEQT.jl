
const INITIAL_READ_SIZE = 1024
const MAX_READ_SIZE     = 16*1024*1024

@assert MAX_READ_SIZE%INITIAL_READ_SIZE==0 && ispow2(div(MAX_READ_SIZE,INITIAL_READ_SIZE)) "MAX_READ_SIZE should be a power of two multiple of INITIAL_READ_SIZE"


struct Tag
	t::UInt16 # two-character representation
end
Tag(a, b) = Tag(UInt16(a)&0xff | (UInt16(b)<<8))
function Tag(s::String)
	@assert length(s)==2 "Length of tag string must be 2"
	Tag(s[1],s[2])
	#Tag(UInt16(s[1])&0xff | (UInt16(s[2])&0xff)<<8)
end

# allow writing tag"NH"
macro tag_str(p)
  Tag(p)
end

print(io::IO,tag::Tag) = print(io,Char(tag.t&0xff),Char(tag.t>>8))
show(io::IO,tag::Tag) = print(io,tag)



mutable struct BamRead
	file::BamFile # file containing the read

	# fields of known size
	block_size::Int32
	refID::Int32
	pos::Int32
	bin_mq_nl::UInt32
	flag_nc::UInt32
	l_seq::Int32
	next_refID::Int32
	next_pos::Int32
	tlen::Int32

	# remaining data of unknown size
	buf::Array{UInt8,1}

	# information to make indexing into buf quicker
	#cigar_start::Int
	seq_start::Int
	qual_start::Int
	#tags_start::Int
	

	# tags_parsed::Bool
	# tagDict::Dict{Tag,Int} # TAG->offset to start of tag (type byte) in buf
	tagsParsePos::Int # byte offset of next tag to parse
	tagInfo::Array{Tuple{Tag,Int},1}
end

#BamRead() = BamRead(0,0,0,0,0,0,0,0,0,Array{UInt8,1}(INITIAL_READ_SIZE),-1,-1)
#BamRead() = BamRead(0,0,0,0,0,0,0,0,0,Array{UInt8,1}(INITIAL_READ_SIZE),-1,-1,false,Dict{UInt16,Int}())
# BamRead(f::BamFile) = BamRead(f,0,0,0,0,0,0,0,0,0,Array{UInt8,1}(INITIAL_READ_SIZE),-1,-1,false,Dict{UInt16,Int}())
BamRead(f::BamFile) = BamRead(f,0,0,0,0,0,0,0,0,0,zeros(UInt8,INITIAL_READ_SIZE),-1,-1,0,Tuple{Tag,Int}[])


block_size(r::BamRead) = Int(r.block_size)


refID(r::BamRead) = Int(r.refID)
pos(r::BamRead) = Int(r.pos)+1
bin_mq_nl(r::BamRead) = UInt(r.bin_mq_nl)
flag_nc(r::BamRead) = UInt(r.flag_nc)
l_seq(r::BamRead) = Int(r.l_seq)
next_refID(r::BamRead) = Int(r.next_refID)
next_pos(r::BamRead) = Int(r.next_pos)+1
tlen(r::BamRead) = Int(r.tlen)


# access to things encoded in the fields
bin(r::BamRead)         = Int(r.bin_mq_nl>>16)
mapq(r::BamRead)        = Int((r.bin_mq_nl>>8)&0x000000ff)
l_read_name(r::BamRead) = Int(r.bin_mq_nl&0x000000ff)
flag(r::BamRead)        = Int(r.flag_nc>>16)
n_cigar_op(r::BamRead)  = Int(r.flag_nc&0x0000ffff)

flagset(r::BamRead,f::Int64) = (flag(r)&f) != 0


# access to things using BamFile info
ref_name(r::BamRead)      = sequence_name(r.file,refID(r))
next_ref_name(r::BamRead) = sequence_name(r.file,next_refID(r))



getuint32(r::BamRead, byteOffs::Int) = (UInt32(r.buf[byteOffs])) | (UInt32(r.buf[byteOffs+1])<<8) | (UInt32(r.buf[byteOffs+2])<<16) | (UInt32(r.buf[byteOffs+3])<<24)
getuint16(r::BamRead, byteOffs::Int) = (UInt16(r.buf[byteOffs])) | (UInt16(r.buf[byteOffs+1])<<8)
getuint8( r::BamRead, byteOffs::Int) = r.buf[byteOffs]

getint32(r::BamRead, byteOffs::Int) = reinterpret(Int32, getuint32(r,byteOffs))
getint16(r::BamRead, byteOffs::Int) = reinterpret(Int16, getuint16(r,byteOffs))
getint8( r::BamRead, byteOffs::Int) = reinterpret(Int8, getuint8(r,byteOffs))


function compute_offsets!(r::BamRead)
	#r.cigar_offs = l_read_name(r)
	#r.seq_offs   = r.cigar_offs + 4*n_cigar_op(r)
	r.seq_start   = l_read_name(r)+1 + 4*n_cigar_op(r)
	r.qual_start  = r.seq_start + ((l_seq(r)+1)>>1)
	#r.tags_start  = r.qual_start + l_seq(r)
	# r.tags_parsed = false
	# isempty(r.tagDict) || (r.tagDict = Dict{Tag,Int}()) # TODO: avoid creating a new dict for every read, this causes too many allocations.

	r.tagsParsePos = r.qual_start + l_seq(r)
	resize!(r.tagInfo,0)

	return
end


# overwrite read with data from stream
#function nextread!(io::IO,r::BamRead)
function nextread!(io,r::BamRead)
	r.block_size = read(io, Int32)
	r.refID      = read(io, Int32)
	r.pos        = read(io, Int32)
	r.bin_mq_nl  = read(io, UInt32)
	r.flag_nc    = read(io, UInt32)
	r.l_seq      = read(io, Int32)
	r.next_refID = read(io, Int32)
	r.next_pos   = read(io, Int32)
	r.tlen       = read(io, Int32)

	reqsize = Int(r.block_size)-8*4
	#reqsize > length(r.buf) && error("Read buffer too small!") # we could expand it, but if the file is corrupt, we could allocate huge things here...
	@assert reqsize>=0 "Corrupt Bam Read, block_size too small"

	# expand read buffer within reasonable limits
	# println("pointer: ", pointer(r.buf), ", size: ", length(r.buf))
	if reqsize>length(r.buf)
		sz = nextpow2(reqsize)
		@assert sz<=MAX_READ_SIZE "Read does not fit in maximum read buffer size ($MAX_READ_SIZE bytes)."
		r.buf = Array{UInt8,1}(sz)
	end

	# println("reqsize: ", reqsize)

	n = readbytes!(io,r.buf,reqsize)
	n!=reqsize && error("Unexpected end of stream when loading read.")
	compute_offsets!(r)
	r
end





# CIGAR and Seq needs to be (immutable) types rather than simple arrays so that we can interpret them for the user

# TODO: should CIGAR be a subtype of AbstractArray???
struct CIGAR{T<:AbstractArray{UInt8}}
	# TODO: Probably make the array type a type parameter. 
	# This way, it also makes sure we don't need to wrap the array into a subarray when it's not necessary.
	# buf::SubArray{UInt8,1,Array{UInt8,1},Tuple{UnitRange{Int}},true}
	buf::T
end

# allow user to construct CIGAR from CIGAR string
function CIGAR(str::AbstractString)
	# len = 0 # compute length of CIGAR
	# for c in str
	# 	len += isalpha(c)
	# end

	# compute cigar length and preallocate array
	cigarLength = length(str) - mapreduce(isnumeric,+,str)
	buf = zeros(UInt32,cigarLength)

	# now interpret cigar string
	len = 0
	ops = Dict([('M',CIGAR_M), ('I',CIGAR_I), ('D',CIGAR_D), ('N',CIGAR_N), ('S',CIGAR_S), ('H',CIGAR_H), ('P',CIGAR_P), ('=',CIGAR_Eq), ('X',CIGAR_X)])

	i = 1
	for c in str
		if isnumeric(c)
			len = len*10 + (c-'0')
			continue
		end
		@assert len>0 "Malformed CIGAR string."
		buf[i] = (len<<4) | ops[c]
		len = 0
		i = i+1
	end
	CIGAR(view(reinterpret(UInt8,buf),1:cigarLength*4))
end


struct SeqElement
	e::Int
end

# TODO: should Seq be a subtype of AbstractArray???
struct Seq <: AbstractArray{SeqElement,1}
	#buf::SubArray{UInt8,1,Array{UInt8,1},Tuple{UnitRange{Int}},1}
	buf::Array{UInt8,1} # the entire buffer with the read (same as BamRead.buf)
	#start::Int          # offset into read buffer - TODO: change to half-bytes to allow substrings???
	start::Int         # offset in half-bytes into read buffer
	len::Int            # length in half-bytes!
end


# access to variable length fields
# read_name(r::BamRead) = SubString(String(r.buf), 1, l_read_name(r)-1)
read_name(r::BamRead) = String(r.buf[1:l_read_name(r)-1])
cigar(r::BamRead)     = (s=l_read_name(r)+1; CIGAR(view(r.buf, s:s+n_cigar_op(r)*4-1)))
seq(r::BamRead)       = Seq(r.buf, (r.seq_start-1)*2, l_seq(r))
qual(r::BamRead)      = view(r.buf, r.qual_start:r.qual_start+l_seq(r)-1)



# TODO: make sure tags are properly copied if we allow the read to not start at the beginning of buf
function copy(r::BamRead)
	BamRead(r.file, r.block_size, r.refID, r.pos, r.bin_mq_nl, r.flag_nc, 
			r.l_seq, r.next_refID, r.next_pos, r.tlen, 
			copy(view(r.buf,1:r.block_size-8*4)), 
			r.seq_start, r.qual_start, r.tagsParsePos, copy(r.tagInfo))
end


copy(c::CIGAR) = CIGAR(view(copy(c.buf),1:length(c.buf))) # copy to array and make SubArray covering the whole
#copy(s::Seq)   = Seq(copy(view(s.buffer,s.start:s.start+((s.len+1)>>1))), 1, s.len)
copy(s::Seq)   = Seq(copy(view(s.buf,s.start>>1+1:(s.start+s.len+1)>>1)), mod(s.start,2), s.len)





struct CIGARElement
	e::UInt32
end
CIGARElement(c::CIGAR,i::Int) = CIGARElement((UInt32(c.buf[i])) | (UInt32(c.buf[i+1])<<8) | (UInt32(c.buf[i+2])<<16) | (UInt32(c.buf[i+3])<<24))

op(e::CIGARElement) = Int(e.e&0xf)
const op_chars = Char['M','I','D','N','S','H','P','=','X']
op_char(e::CIGARElement) = op_chars[op(e)+1]
len(e::CIGARElement) = Int(e.e>>4)


print(io::IO,e::CIGARElement) = print(io,len(e),op_char(e))



# CIGAR iterators
# Iterators
# start(::CIGAR) = 1
# next(c::CIGAR, i::Int) = (CIGARElement(c,i), i+4)
# done(c::CIGAR, i::Int) = i > length(c.buf)

function iterate(c::CIGAR, i::Int=1)
	i>length(c.buf) && return nothing
	(CIGARElement(c,i), i+4)
end
eltype(::Type{T}) where {T<:CIGAR} = CIGARElement
length(c::CIGAR) = div(length(c.buf),4)



function alignedlength(cig::CIGAR)
	l = 0
	for c in cig
		o = op(c)
		if o==CIGAR_M || o==CIGAR_X || o==CIGAR_Eq || o==CIGAR_D || o==CIGAR_N
			l += len(c)
		end
	end
	l
end


endpos(r::BamRead) = pos(r) + alignedlength(cigar(r)) - 1




function SeqElement(s::Seq, i::Int)
	b = Int(s.buf[(i+1)>>1]) # get the correct byte
	b = i&1==0 ? b&0x0f : b>>4 # select low or high nybble depending on index
	SeqElement(b)
end

const base_indices = [0,1,2,0,3,0,0,0,4,0,0,0,0,0,0,0]
const base_chars = UInt8['=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'];
base_mask(e::SeqElement) = e.e
base_index(e::SeqElement) = base_indices[e.e+1]
base_char(e::SeqElement) = base_chars[e.e+1] # TODO: should this be renamed, return Char or change in some other way?


# Iterators
# start(::Seq) = 1
# next(s::Seq, i::Int) = (SeqElement(s,i+s.start), i+1)
# done(s::Seq, i::Int) = i > s.len
function iterate(s::Seq, i::Int=1)
	i > s.len && return nothing
	(SeqElement(s,i+s.start), i+1)
end
eltype(::Type{Seq}) = SeqElement



size(s::Seq) = (s.len,)
length(s::Seq) = s.len
# linearindexing(::Type{Seq}) = Base.LinearFast()
Base.IndexStyle(::Type{Seq}) = Base.IndexStyle(Type{Array})


# TODO: is this the proper way to do bounds checking???
checkbounds(s::Seq, i::Int) = 1<=i<=s.len || throw(BoundsError("$(s.len)-element Seq",i))
checkbounds(::Type{Bool}, s::Seq, i::Int) = 1<=i<=s.len

# random access 
#getindex(s::Seq, i::Int) = (checkbounds(s,i); SeqElement(s,i+(s.start-1)*2))
getindex(s::Seq, i::Int) = (checkbounds(s,i); SeqElement(s,i+s.start))


# construct view into seq
#view(s::Seq, rng::UnitRange{Int}) = (checkbounds(s,rng[1]); checkbounds(s,rng[end]); Seq(s.buf,s.start+rng[1]-1,length(rng)))
function view(s::Seq, rng::UnitRange{Int}) 
	length(rng)==0 && return Seq(s.buf,0,0) # anything with length 0 will do
	checkbounds(s,rng[1])
	checkbounds(s,rng[end])
	Seq(s.buf,s.start+rng[1]-1,length(rng))
end




function print(io::IO,cig::CIGAR)
	for c in cig
		print(io,c)
	end
end
show(io::IO,cig::CIGAR) = print(io,cig)

# TODO: implement this in a better (still typestable) way.
function string(cig::CIGAR)
	io = IOBuffer()
	print(io,cig)
	String(take!(io)) # for type stability
end
cigar_str(r::BamRead) = string(cigar(r)) # convenience access



function print(io::IO,seq::Seq)
	for (i,n) in enumerate(seq)
		print(io,Char(base_char(n)))
	end
end
show(io::IO,seq::Seq) = print(io,seq)
function string(seq::Seq)
	str = zeros(UInt8,seq.len) # allocate buffer
	for (i,n) in enumerate(seq)
		str[i] = base_char(n)
	end
	String(str)
end
seq_str(r::BamRead) = string(seq(r)) # convenience access




# tags - array implementation


read_tag_name(r::BamRead, offs::Int) = Tag(r.buf[offs], r.buf[offs+1])
read_tag_type(r::BamRead, offs::Int) = Char(r.buf[offs])

# sizes for fixed-size tags
const tagTypeSizes = Dict{Char,Int}('A'=>1,'c'=>1,'C'=>1,'s'=>2,'S'=>2,'i'=>4,'I'=>4,'f'=>4)


function parsenexttag!(r::BamRead)
	tag = read_tag_name(r,r.tagsParsePos) # two-character tag
	r.tagsParsePos += 2

	# store offset in array
	push!(r.tagInfo, (tag, r.tagsParsePos))

	# now skip rest of tag
	tagType = read_tag_type(r,r.tagsParsePos)
	r.tagsParsePos += 1

	if tagType=='Z'
		# iterate until 0 character
		while r.buf[r.tagsParsePos] != 0x0
			r.tagsParsePos += 1
		end
		r.tagsParsePos += 1
	elseif tagType=='B'
		error("Tags with array type (B) are not yet supported.")
	elseif tagType=='H' # is this tag possible in BAM???
		error("Tags with type H are not yet supported.")
	else
		# fixed-size tag
		r.tagsParsePos += tagTypeSizes[tagType]
	end
	tag
end


# find tag and parse offset until it is found
function findtag!(r::BamRead, tag::Tag)

	# first check in already parsed tags
	for (currTag, offs) in r.tagInfo
		currTag == tag && return offs
	end


	# and if it wasn't found, continue looking in unparsed part
	while r.tagsParsePos < r.block_size-8*4 # used part of buf
		currTag = parsenexttag!(r)
		currTag == tag && return r.tagInfo[end][2]
	end

	-1
end



hastag(r::BamRead, tag::Tag) = findtag!(r,tag) != -1
hastag(r::BamRead, tag::String) = hastag(r,Tag(tag))

# TODO: support default value if tag doesn't exist?
function tagint(r::BamRead, tag::Tag)
	offs = findtag!(r,tag)
	tagType = read_tag_type(r,offs)
	offs += 1

	#println("tagtype: $tagType")

	# get value from byte buffer
	if tagType=='i'
		v = Int(getint32(r,offs))
	elseif tagType=='I'
		v = Int(getuint32(r,offs))
	elseif tagType=='s'
		v = Int(getint16(r,offs))
	elseif tagType=='S'
		v = Int(getuint16(r,offs))
	elseif tagType=='c'
		v = Int(getint8(r,offs))
	elseif tagType=='C'
		v = Int(getuint8(r,offs))
	else
		error("Cannot convert tag type \"$tagType\" to Int")
	end
	v
end
tagint(r::BamRead, tag::String) = tagint(r,Tag(tag))

# TODO: support default value if tag doesn't exist?
function tagchar(r::BamRead, tag::Tag)
	offs = findtag!(r,tag)
	tagType = read_tag_type(r,offs)
	offs += 1

	@assert tagType=='A' "Cannot convert tag type \"$tagType\" to Char"
	Char(getuint8(r,offs))
end
tagchar(r::BamRead, tag::String) = tagchar(r,Tag(tag))



# TODO: support default value if tag doesn't exist?
function tagfloat32(r::BamRead, tag::Tag)
	offs = findtag!(r,tag)
	tagType = read_tag_type(r,offs)
	offs += 1

	@assert tagType=='f' "Cannot convert tag type \"$tagType\" to Float64"
	reinterpret(Float32,getuint32(r,offs))
end
tagfloat32(r::BamRead, tag::String) = tagfloat(r,Tag(tag))

tagfloat(r::BamRead, tag) = Float(tagfloat(r,tag))


# TODO: support default value if tag doesn't exist?
function tagstring(r::BamRead, tag::Tag)
	offs = findtag!(r,tag)
	tagType = read_tag_type(r,offs)
	offs += 1

	@assert tagType=='Z' "Cannot convert tag type \"$tagType\" to string"

	start = offs

	# find end of string - TODO: no need to do this twice, store end pos in dict somehow?
	while r.buf[offs] != 0x0
		offs += 1
	end
	offs -= 1 # do not include NULL terminator
	
	# SubString(String(r.buf),start,offs)
	String(r.buf[start:offs])
end
tagstring(r::BamRead, tag::String) = tagstring(r, Tag(tag))


function printtag(io::IO, r::BamRead, tag::Tag, offs::Int)
	tagType = read_tag_type(r,offs)
	if tagType=='i' || tagType=='s' || tagType=='c' || tagType=='I' || tagType=='S' || tagType=='C'
		print(io, '\t', tag, ":i:", tagint(r,tag))
	elseif tagType=='A'
		print(io, '\t', tag, ":A:", tagchar(r,tag))
	elseif tagType=='f'
		print(io, '\t', tag, ":f:", tagfloat(r,tag)) # TODO: check exactly how floats are normally printed to SAM
	elseif tagType=='Z'
		print(io, '\t', tag, ":Z:", tagstring(r,tag))
	else
		error("Tags with array type (", tagType, ") are not currently supported.")
	end	
end

function printtags(io::IO, r::BamRead)
	# first loop over already parsed tags
	for (tag, offs) in r.tagInfo
		printtag(io, r, tag, offs)
	end

	# and continue with previously unparsed tags
	while r.tagsParsePos < r.block_size-8*4 # used part of buf
		tag = parsenexttag!(r)
		printtag(io, r, tag, r.tagInfo[end][2])
	end
end



# # tags - dict implementation

# #read_tag_name(r::BamRead, offs::Int) = UInt16(r.buf[offs]) | UInt16(r.buf[offs+1])<<8
# read_tag_name(r::BamRead, offs::Int) = Tag(r.buf[offs], r.buf[offs+1])
# read_tag_type(r::BamRead, offs::Int) = Char(r.buf[offs])

# # sizes for fixed-size tags
# const tagTypeSizes = Dict{Char,Int}('A'=>1,'c'=>1,'C'=>1,'s'=>2,'S'=>2,'i'=>4,'I'=>4,'f'=>4)

# function parsetags(r::BamRead)
# 	offs = r.qual_start + l_seq(r)

# 	while offs < r.block_size-8*4 # used part of buf
# 		tag = read_tag_name(r,offs) # two-character tag
# 		offs += 2

# 		# store offset in dict
# 		r.tagDict[tag] = offs 

# 		# now skip rest of tag
# 		tagType = read_tag_type(r,offs)
# 		offs += 1

# 		if tagType=='Z'
# 			#error("Tags with string type are not yet supported.")
# 			# println(SubString(String(r.buf),offs,offs+10-1))
# 			# println(r.buf[offs+10])

# 			# iterate until 0 character
# 			while r.buf[offs] != 0x0
# 				offs += 1
# 			end
# 			offs += 1
# 		elseif tagType=='B'
# 			error("Tags with array type (B) are not yet supported.")
# 		elseif tagType=='H' # is this tag possible in BAM???
# 			error("Tags with type H are not yet supported.")
# 		else
# 			# fixed-size tag
# 			offs += tagTypeSizes[tagType]
# 		end
# 	end


# 	r.tags_parsed = true
# 	return
# end


# function hastag(r::BamRead, tag::Tag)
# 	r.tags_parsed || parsetags(r)
# 	#tagname2uint16(tag) ∈ keys(r.tagDict)
# 	haskey(r.tagDict,tag)
# end
# hastag(r::BamRead, tag::String) = hastag(r,Tag(tag))

# # TODO: support default value if tag doesn't exist?
# function tagint(r::BamRead, tag::Tag)
# 	r.tags_parsed || parsetags(r)
# 	offs = r.tagDict[tag]
# 	tagType = read_tag_type(r,offs)
# 	offs += 1

# 	#println("tagtype: $tagType")

# 	# get value from byte buffer
# 	if tagType=='i'
# 		v = Int(getint32(r,offs))
# 	elseif tagType=='I'
# 		v = Int(getuint32(r,offs))
# 	elseif tagType=='s'
# 		v = Int(getint16(r,offs))
# 	elseif tagType=='S'
# 		v = Int(getuint16(r,offs))
# 	elseif tagType=='c'
# 		v = Int(getint8(r,offs))
# 	elseif tagType=='C'
# 		v = Int(getuint8(r,offs))
# 	else
# 		error("Cannot convert tag type \"$tagType\" to Int")
# 	end
# 	v
# end
# tagint(r::BamRead, tag::String) = tagint(r,Tag(tag))

# # TODO: support default value if tag doesn't exist?
# function tagchar(r::BamRead, tag::Tag)
# 	r.tags_parsed || parsetags(r)
# 	offs = r.tagDict[tag]
# 	tagType = read_tag_type(r,offs)
# 	offs += 1

# 	@assert tagType=='A' "Cannot convert tag type \"$tagType\" to Char"

# 	Char(getuint8(r,offs))
# end
# tagchar(r::BamRead, tag::String) = tagchar(r,Tag(tag))



# # TODO: support default value if tag doesn't exist?
# function tagstring(r::BamRead, tag::Tag)
# 	r.tags_parsed || parsetags(r)
# 	offs = r.tagDict[tag]
# 	tagType = read_tag_type(r,offs)
# 	offs += 1

# 	@assert tagType=='Z' "Cannot convert tag type \"$tagType\" to string"

# 	start = offs

# 	# find end of string - TODO: no need to do this twice, store end pos in dict somehow?
# 	while r.buf[offs] != 0x0
# 		offs += 1
# 	end
# 	offs -= 1 # do not include NULL terminator
	
# 	SubString(String(r.buf),start,offs)
# end
# tagstring(r::BamRead, tag::String) = tagstring(r, Tag(tag))





# create SAM string representation of BamRead
function print(io::IO, r::BamRead)
	# QNAME
	# FLAG
	# RNAME
	# POS
	# MAPQ
	# CIGAR
	# RNEXT
	# PNEXT
	# TLEN
	# SEQ
	# QUAL
	print(io, read_name(r), '\t')
	print(io, flag(r), '\t')
	print(io, ref_name(r), '\t')
	print(io, pos(r), '\t')
	print(io, mapq(r), '\t')
	cig = string(cigar(r))
	print(io, length(cig)==0 ? '*' : cig, '\t')
	#print(io, next_ref_name(r), '\t')
	print(io, next_refID(r)==refID(r) ? "=" : next_ref_name(r), '\t')
	print(io, next_pos(r), '\t')
	print(io, tlen(r), '\t')
	print(io, seq(r), '\t')
	print(io, String(qual(r).+0x21)) # NO tab, will be added by tag printing if necessary
	
	printtags(io, r)
end
show(io::IO, r::BamRead) = print(io,r)


function sam(r::BamRead)
	io = IOBuffer()
	print(io,r)
	String(take!(io))
end

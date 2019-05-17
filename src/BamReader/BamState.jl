import Base: iterate, eltype

# very simple implementation of iterating through all the reads of a BAM file

struct BamState
	io::IO
	zio::IO
	read::BamRead
end

BamState(f::BamFile, io::IO, zio::IO) = BamState(io,zio,BamRead(f))
function BamState(f::BamFile, io::IO) 
	zio = ZlibInflateInputStream(io,reset_on_end=true)
	skip_bam_header(zio) # just skip header for now
	BamState(f, io, zio)
end
BamState(f::BamFile) = BamState(f,open(f.filename))



# Base.start(f::BamFile) = BamState(f)
# Base.next(::BamFile, s::BamState) = (nextread!(s.zio, s.read), s)
# Base.done(::BamFile, s::BamState) = eof(s.zio)


function iterate(f::BamFile, s::BamState=BamState(f))
	eof(s.zio) && return nothing
	(nextread!(s.zio, s.read), s)
end
Base.eltype(::Type{BamFile}) = BamRead

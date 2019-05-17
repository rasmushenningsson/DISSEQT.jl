

function read_bam_header(io::IO)
	magic = read(io,UInt32)
	magic != 0x014d4142 && error("BAM magic string not found")

	l_text = read(io,Int32)
	# text = String(read(io,UInt8,l_text))
	text = String(read!(io,Vector{UInt8}(undef,l_text)))
	n_ref = read(io,Int32)

	refs = Vector{Tuple{String,Int}}(undef,n_ref)

	for i=1:n_ref
		l_name = read(io,Int32)
		# name = read(io,UInt8,l_name)
		name = read!(io,Vector{UInt8}(undef,l_name))
		l_ref = read(io,Int32)
		refs[i] = (String(name[1:end-1]),l_ref)
	end
	(text,refs)
end
function read_bam_header(filename)
	io = open(filename)	
	zio = ZlibInflateInputStream(io,reset_on_end=true)
	(text,refs) = read_bam_header(zio)
	#close(zio)
	close(io)
	(text,refs)
end

function skip_bam_header(io::IO)
	magic = read(io,UInt32)
	magic != 0x014d4142 && error("BAM magic string not found")

	l_text = read(io,Int32)
	text = read!(io,Vector{UInt8}(undef,l_text))
	n_ref = read(io,Int32)

	for i=1:n_ref
		l_name = read(io,Int32)
		name = read!(io,Vector{UInt8}(undef,l_name))
		l_ref = read(io,Int32)
	end
end


mutable struct BamFile
	filename::AbstractString

	# header
	text::String
	refs::Vector{Tuple{String,Int}}

	# index
	index::BamIndex
end
BamFile() = BamFile("","",[],BamIndex())
BamFile(filename) = BamFile(filename, read_bam_header(filename)..., tryloadindex(filename))


filename(f::BamFile) = f.filename
header(f::BamFile) = f.text
sequences(f::BamFile) = f.refs
sequence_name(f::BamFile, refID::Integer) = refID<0 ? "*" : f.refs[refID+1][1]

refID(f::BamFile, name::AbstractString) = findfirst(x->x[1]==name, sequences(f))-1

hasindex(f::BamFile) = !isempty(f.index)
index(f::BamFile) = f.index

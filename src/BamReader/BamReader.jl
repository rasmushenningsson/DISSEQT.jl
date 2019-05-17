module BamReader

using DataStructures
using Libz
using BufferedStreams

import Base: show, print, seek, +, isless, isempty, copy, iterate, eltype, size, length,
             getindex, checkbounds, string, view, position, eof, read, readbytes!



export 
	BamFile,
	BamRead,
	CIGAR,
	#CIGARState,
	CIGARElement,
	Seq,
	#SeqState,
	SeqElement,
	Quality,
	QualityState,
	QualityElement,
	filename,
	header,
	sequences,
	sequence_name,

	refID,
	pos,
	bin_mq_nl,
	flag_nc,
	l_seq,
	next_refID,
	next_pos,
	tlen,
	
	# bin, # TODO: find a name that doesn't clash with Base.bin
	mapq,
	l_read_name,
	flag,
	n_cigar_op,
	flagset,

	ref_name,
	next_ref_name,
	
	read_name,
	cigar,
	seq,
	qual,

	op,
	op_char,
	len,
	cigar_str,
	alignedlength,
	endpos,
	base_mask,
	base_index,
	base_char,
	seq_str,
	# quality_array,
	@tag_str,
	hastag,
	tagint,
	tagchar,
	tagfloat,
	tagstring,

	sam, 

	BAM_FPAIRED,
	BAM_FPROPER_PAIR,
	BAM_FUNMAP,
	BAM_FMUNMAP,
	BAM_FREVERSE,
	BAM_FMREVERSE,
	BAM_FREAD1,
	BAM_FREAD2,
	BAM_FSECONDARY,
	BAM_FQCFAIL,
	BAM_FDUP,
	BAM_FSUPPLEMENTARY,
	CIGAR_M,
	CIGAR_I,
	CIGAR_D,
	CIGAR_N,
	CIGAR_S,
	CIGAR_H,
	CIGAR_P,
	CIGAR_Eq,
	CIGAR_X,

	readpos,
	refpos,

	BamIndex,
	reg2bin,
	reg2bins,
	reg2binranges,

	region


include("constants.jl")
include("BamIndex.jl")
include("BamFile.jl")
include("BamRead.jl")
include("BamState.jl")
include("BGZFBlockStream.jl")
include("BamRegion.jl")
include("readiterators.jl")

end
using Test, DISSEQT.BamReader

# TODO: add Seq tests



# usage
# @test @typestable(expr) == Int
# @test @typestable(expr) == subtypeof(Integer)
macro typestable(expr)
	esc(:(typeof(@inferred($expr))))
end


# a little hack to be able to use --> subtypeof(ParentType)
import Base: ==, show
struct SubTypeOf
	t::Union{DataType,UnionAll}
end
==(x::Union{DataType,UnionAll},y::SubTypeOf) = x <: y.t
==(x::SubTypeOf,y::Union{DataType,UnionAll}) = y==x
show(io::IO, s::SubTypeOf) = print(io,"<:",s.t)

subtypeof(t::Union{DataType,UnionAll}) = SubTypeOf(t)



@testset "CIGAR" begin
	strs = ["1M", "10S268435455M21S", "20S10M1I9M2D10M3X7=10M4P24M10H"]
	cigArrays = Array[[CIGARElement(1<<4 | CIGAR_M)],
	                 [CIGARElement(10<<4 | CIGAR_S), CIGARElement(268435455<<4 | CIGAR_M), CIGARElement(21<<4 | CIGAR_S)],
				     [CIGARElement(20<<4 | CIGAR_S), CIGARElement(10<<4 | CIGAR_M), CIGARElement(1<<4 | CIGAR_I), CIGARElement(9<<4 | CIGAR_M), CIGARElement(2<<4 | CIGAR_D), CIGARElement(10<<4 | CIGAR_M), CIGARElement(3<<4 | CIGAR_X), CIGARElement(7<<4 | CIGAR_Eq), CIGARElement(10<<4 | CIGAR_M), CIGARElement(4<<4 | CIGAR_P), CIGARElement(24<<4 | CIGAR_M), CIGARElement(10<<4 | CIGAR_H)]]
	for (str,cigArray) in zip(strs,cigArrays)
		cig = CIGAR(str)
		@test string(cig) == str
		@test collect(cig) == cigArray
	end
end


function verifymyread(r::BamRead)
	rname = "TheReadName";

	@test refID(r) == 0
	@test @typestable(refID(r)) == Int

	@test pos(r) == 1234
	@test @typestable(pos(r)) == Int

	@test bin_mq_nl(r) == UInt((Int(0x1249)<<16)|(14<<8)|(length(rname)+1))
	@test @typestable(bin_mq_nl(r)) == UInt

	@test flag_nc(r) == UInt(0x600000c)
	@test @typestable(flag_nc(r)) == UInt

	@test l_seq(r) == 94
	@test @typestable(l_seq(r)) == Int

	@test next_refID(r) == -1
	@test @typestable(next_refID(r)) == Int

	@test next_pos(r) == 0
	@test @typestable(next_pos(r)) == Int

	@test tlen(r) == 0
	@test @typestable(tlen(r)) == Int

	@test BamReader.bin(r) == 0x1249
	@test @typestable(BamReader.bin(r)) == Int

	@test mapq(r) == UInt8(14)
	@test @typestable(mapq(r)) == Int

	@test l_read_name(r) == UInt8(length(rname)+1)
	@test @typestable(l_read_name(r)) == Int

	@test flag(r) == 0x600
	@test @typestable(flag(r)) == Int

	@test n_cigar_op(r) == 12
	@test @typestable(n_cigar_op(r)) == Int

	@test read_name(r) == rname
	@test @typestable(read_name(r)) == subtypeof(AbstractString)

	@test ref_name(r) == "TheChrName"
	@test @typestable(ref_name(r)) == subtypeof(AbstractString)

	cig = "20S10M1I9M2D10M3X7=10M4P24M10H"
	@test string(cigar(r)) == cig
	@test @typestable(string(cigar(r))) == subtypeof(AbstractString)

	@test string(copy(cigar(r))) == cig
	@test @typestable(string(copy(cigar(r)))) == subtypeof(AbstractString)

	@test cigar_str(r) == cig
	@test @typestable(cigar_str(r)) == subtypeof(AbstractString)


	cigArray = [CIGARElement(20<<4 | CIGAR_S), CIGARElement(10<<4 | CIGAR_M), CIGARElement(1<<4 | CIGAR_I), CIGARElement(9<<4 | CIGAR_M), CIGARElement(2<<4 | CIGAR_D), CIGARElement(10<<4 | CIGAR_M), CIGARElement(3<<4 | CIGAR_X), CIGARElement(7<<4 | CIGAR_Eq), CIGARElement(10<<4 | CIGAR_M), CIGARElement(4<<4 | CIGAR_P), CIGARElement(24<<4 | CIGAR_M), CIGARElement(10<<4 | CIGAR_H)]
	@test collect(cigar(r)) == cigArray
	@test @typestable(collect(cigar(r))) == Array{CIGARElement,1}


	s = "ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTTAAAAACCCCCGGGGGTTTTTAAAAAACCCCCCGGGGGGTTTTTTGTCACATTAT"
	@test string(seq(r)) == s
	@test @typestable(string(seq(r))) == subtypeof(AbstractString)

	@test string(copy(seq(r))) == s
	@test @typestable(string(copy(seq(r)))) == subtypeof(AbstractString)

	@test seq_str(r) == s
	@test @typestable(seq_str(r)) == subtypeof(AbstractString)


	@test join(map(x->Char(base_char(x)), seq(r))) == s

	rngs = UnitRange{Int}[ 1:length(s), 2:length(s)-1, 1:length(s)-1, 2:length(s)-2, 3:3, 4:4, 1:2, 2:3, 4:87, 1:0, 2:1, length(s):length(s)-1 ]
	for rng in rngs
		@test string(view(seq(r),rng)) == s[rng]
		@test @typestable(string(view(seq(r),rng))) == subtypeof(AbstractString)
	end


	@test qual(r) == collect(UInt8('!'):UInt8('~')).-33
	@test @typestable(qual(r)) == subtypeof(AbstractArray)

	@test hastag(r,"CC") == true
	@test hastag(r,"AA") == false
	@test tagint(r,"NH") == 1
	@test tagint(r,"CP") == 6789
	@test tagstring(r,"CC") == "AnotherChr"

	@test @typestable(hastag(r,"CC")) == Bool
	@test @typestable(tagint(r,"NH")) == Int
	@test @typestable(tagstring(r,"CC")) == subtypeof(AbstractString)

	@test hastag(r,tag"CC") == true
	@test hastag(r,tag"AA") == false
	@test tagint(r,tag"NH") == 1
	@test tagint(r,tag"CP") == 6789
	@test tagstring(r,tag"CC") == "AnotherChr"

	@test @typestable(hastag(r,tag"CC")) == Bool
	@test @typestable(tagint(r,tag"NH")) == Int
	@test @typestable(tagstring(r,tag"CC")) == subtypeof(AbstractString)


	# TODO: test _all_ tag types


	# compare printout with SAM file
	samfromfile = read("BamReader/data/singleread.sam", String)
	samfromfile = split(samfromfile,"\n")[3] # hardcoded that there are two header rows
	@test sam(r) == samfromfile
end



@testset "BamRead" begin
	buf = UInt8[0,1,0,0,0,0,0,0,209,4,0,0,12,14,73,18,12,0,0,6,94,0,0,0,255,255,255,255,255,255,255,255,0,0,0,0,84,104,101,82,101,97,100,78,97,109,101,0,68,1,0,0,160,0,0,0,17,0,0,0,144,0,0,0,34,0,0,0,160,0,0,0,56,0,0,0,119,0,0,0,160,0,0,0,70,0,0,0,128,1,0,0,165,0,0,0,18,72,17,34,68,136,17,18,34,68,72,136,17,17,34,34,68,68,136,136,17,17,18,34,34,68,68,72,136,136,17,17,17,34,34,34,68,68,68,136,136,136,72,33,33,136,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,78,72,67,1,67,67,90,65,110,111,116,104,101,114,67,104,114,0,67,80,83,133,26]

	bamFile = BamFile("dummyname", "@HD\tVN1.3\n@SQ\tSN:TheChrName\tLN:1000000\n", [("TheChrName",1000000)], BamIndex())
	
	# verify header
	@test filename(bamFile) == "dummyname"
	@test header(bamFile) == "@HD\tVN1.3\n@SQ\tSN:TheChrName\tLN:1000000\n"
	@test sequences(bamFile) == [("TheChrName",1000000)]

	rd = BamRead(bamFile)
	BamReader.nextread!(IOBuffer(buf),rd)
	BamReader.compute_offsets!(rd)

	for r in [rd,copy(rd)] # test on both read and copy of read
		verifymyread(r)
	end
end


# test with small BAM
@testset "BAM" begin
	#bamData = UInt8[0,0,0,0,209,4,0,0,12,14,73,18,12,0,71,202,94,0,0,0,255,255,255,255,255,255,255,255,0,0,0,0,84,104,101,82,101,97,100,78,97,109,101,0,68,1,0,0,160,0,0,0,17,0,0,0,144,0,0,0,34,0,0,0,160,0,0,0,56,0,0,0,119,0,0,0,160,0,0,0,70,0,0,0,128,1,0,0,165,0,0,0,18,72,17,34,68,136,17,18,34,68,72,136,17,17,34,34,68,68,136,136,17,17,18,34,34,68,68,72,136,136,17,17,17,34,34,34,68,68,68,136,136,136,72,33,33,136,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,78,72,67,1,67,67,90,65,110,111,116,104,101,114,67,104,114,0,67,80,83,133,26,96,119,5,132,0,0,0,0,192,118,5,132,0,0,0,0,0,120,5,132,0,0,0,0,16,120,5,132,0,0,0,0,48,120,5,132,0,0,0,0,64,120,5,132,0,0,0,0,80,120,5,132,0,0,0,0,240,120,5,132,0,0,0,0,32,121,5,132,0,0,0,0,64,121,5,132,0,0,0,0,80,121,5,132,0,0,0,0,32,123,5,132,0,0,0,0,48,123,5,132,0,0,0,0,48,124,5,132,0,0,0,0,16,123,5,132,0,0,0,0,0,132,6,132,0,0,0,0,32,132,6,132,0,0,0,0,48,134,6,132,0,0,0,0,144,134,6,132,0,0,0,0,160,134,6,132,0,0,0,0,176,134,6,132,0,0,0,0,32,135,6,132,0,0,0,0,160,135,6,132,0,0,0,0,96,137,6,132,0,0,0,0,112,137,6,132,0,0,0,0,48,138,6,132,0,0,0,0,128,138,6,132,0,0,0,0,240,138,6,132,0,0,0,0,0,139,6,132,0,0,0,0,240,139,6,132,0,0,0,0,16,140,6,132,0,0,0,0,112,141,6,132,0,0,0,0,192,141,6,132,0,0,0,0,64,142,6,132,0,0,0,0,48,141,6,132,0,0,0,0,224,154,6,132,0,0,0,0,32,155,6,132,0,0,0,0,96,142,6,132,0,0,0,0,112,142,6,132,0,0,0,0,128,142,6,132,0,0,0,0,160,143,6,132,0,0,0,0,192,143,6,132,0,0,0,0,128,143,6,132,0,0,0,0,192,144,6,132,0,0,0,0,0,145,6,132,0,0,0,0,64,137,6,132,0,0,0,0,16,147,6,132,0,0,0,0,32,147,6,132,0,0,0,0,80,147,6,132,0,0,0,0,112,124,5,132,0,0,0,0,144,124,5,132,0,0,0,0,192,124,5,132,0,0,0,0,32,125,5,132,0,0,0,0,80,125,5,132,0,0,0,0,48,121,5,132,0,0,0,0,96,157,6,132,0,0,0,0,176,157,6,132,0,0,0,0,0,158,6,132,0,0,0,0,48,158,6,132,0,0,0,0,224,158,6,132,0,0,0,0,128,125,5,132,0,0,0,0,176,125,5,132,0,0,0,0,64,126,5,132,0,0,0,0,16,121,5,132,0,0,0,0,240,161,6,132,0,0,0,0,80,162,6,132,0,0,0,0,80,126,5,132,0,0,0,0,144,126,5,132,0,0,0,0,176,126,5,132,0,0,0,0,192,126,5,132,0,0,0,0,224,126,5,132,0,0,0,0,32,127,5,132,0,0,0,0,96,127,5,132,0,0,0,0,16,127,5,132,0,0,0,0,16,128,6,132,0,0,0,0,96,128,6,132,0,0,0,0,64,128,6,132,0,0,0,0,64,130,6,132,0,0,0,0,96,114,5,132,0,0,0,0,80,164,6,132,0,0,0,0,96,164,6,132,0,0,0,0,176,164,6,132,0,0,0,0,192,164,6,132,0,0,0,0,80,165,6,132,0,0,0,0,96,165,6,132,0,0,0,0,208,165,6,132,0,0,0,0,240,165,6,132,0,0,0,0,224,166,6,132,0,0,0,0,112,167,6,132,0,0,0,0,208,168,6,132,0,0,0,0,240,168,6,132,0,0,0,0,160,169,6,132,0,0,0,0,32,170,6,132,0,0,0,0,224,170,6,132,0,0,0,0,0,171,6,132,0,0,0,0,192,175,6,132,0,0,0,0]
	# TODO: work with BAM from memory instead of from file?

	fn = "BamReader/data/singleread.bam"
	bamFile = BamFile(fn)

	# verify header
	@test filename(bamFile) == fn
	@test header(bamFile) == "@HD\tVN1.3\n@SQ\tSN:TheChrName\tLN:1000000\n"
	@test sequences(bamFile) == [("TheChrName",1000000)]


	for (i,rd) in enumerate(bamFile)
		@test i == 1 # there is only one read in this file

		for r in [rd,copy(rd)] # test on both read and copy of read
			verifymyread(r)
		end
	end


end




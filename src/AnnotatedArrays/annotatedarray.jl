
import Base: size, getindex, setindex!, copy, view,
             isequal, ==, !=,# .==, .!=, .<, .<=, .>, .>=, 
             #broadcast,
             <,<=,>,>=,
             # &, |, 
             !, ~,
             cat, vcat, hcat,
             #+, -,
             # transpose, ctranspose,
             convert,
             show, print, display


abstract type AnnotatedIndex end


struct SingletonArray{T,N,S} <: AbstractArray{T,N}
	array::S # where S <: AbstractArray{T,N}
end
SingletonArray(A::AbstractArray{T,N}) where {T,N} = SingletonArray{T,N,typeof(A)}(A)



convert(::Type{Array}, a::SingletonArray) = a.array



isequal(A::SingletonArray,B::SingletonArray) = isequal(A.array,B.array)

==(A::SingletonArray,B::SingletonArray) = A.array==B.array
!=(A::SingletonArray,B::SingletonArray) = A.array!=B.array

broadcast(::typeof(==),A::SingletonArray,B::SingletonArray) = A.array.==B.array
broadcast(::typeof(!=),A::SingletonArray,B::SingletonArray) = A.array.!=B.array
broadcast(::typeof(<),A::SingletonArray,B::SingletonArray) = A.array.>B.array
broadcast(::typeof(<=),A::SingletonArray,B::SingletonArray) = A.array.>=B.array
broadcast(::typeof(>),A::SingletonArray,B::SingletonArray) = A.array.>B.array
broadcast(::typeof(>=),A::SingletonArray,B::SingletonArray) = A.array.>=B.array



size(A::SingletonArray) = size(A.array)
Base.IndexStyle(::Type{SingletonArray}) = Base.IndexStyle(Type{Array})
getindex(A::SingletonArray, i::Int)	= getindex(A.array,i)
setindex!(A::SingletonArray, v, i::Int) = setindex(A.array,v,i)

# TODO: review what happens with empty indices (i.e. [], [false, false, ..., false], etc.)
singletonindex(sz::Int, i::Int) = sz==1 ? 1 : i
singletonindex(sz::Int, r::AbstractRange) = sz==1 ? Colon() : r # TODO: make type stable by creating range of same type with only 1 element
singletonindex(sz::Int, ::Colon) = Colon()

singletonindex(sz::Int, arr::Array{Int,1}) = sz==1 ? [1] : arr

function singletonindex(sz::Int, arr::AbstractArray) 
	@assert isempty(arr) "Can only index with empty array (with element type Any)."
	#arr
	#sz==1 ? [1] : arr # not type stable
	sz==1 ? [1] : Int[] # type stable
end
singletonindex(sz::Int, b::AbstractArray{Bool}) = sz==1 ? trues(1) : b # TODO: make type stable when b is not a BitArray



function getindex(A::SingletonArray, I...)
	# TODO: how do with the number of dimensions???
	sz = size(A,(1:length(I))...) # match number of dimensions of I
	getindex(A.array, map(singletonindex,sz,I)...)
end


function viewimpl(A::SingletonArray, I...)	
	sz = size(A,(1:length(I))...) # match number of dimensions of I
	SingletonArray(view(A.array, map(singletonindex,sz,I)...))
end

view(A::SingletonArray, I::Union{Int64, AbstractVector, Colon}...) = viewimpl(A,I...)


copy(A::SingletonArray) = SingletonArray(copy(A.array))


function cat(dim::Int,A::SingletonArray,A2::SingletonArray...)
	A = (A,A2...)
	SingletonArray(cat(dim, map(x->convert(Array,x),A)...))
end



# transpose(A::SingletonArray) = SingletonArray(transpose(A.array))
# ctranspose(A::SingletonArray) = SingletonArray(ctranspose(A.array))



struct AnnotatedArray{T,N,S} <: AbstractArray{T,N}
	array::S # where S <: AbstractArray{T,N}
	annotations::Dict{Symbol,SingletonArray} # TODO: Use SortedDict for order stability?
end


AnnotatedArray(A::AbstractArray{T,N}) where {T,N} =
	AnnotatedArray{T,N,typeof(A)}(A,Dict{Symbol,SingletonArray}())


convert(::Type{Array}, a::AnnotatedArray) = a.array

annotations(A::AnnotatedArray) = A.annotations
annotationnames(A::AnnotatedArray) = keys(annotations(A))



# the annotation size in each dimension must equal the array size
# or be a singleton dimension
function allowedsize(szArray::Tuple,szAnnotation::Tuple)
	all(i->szAnnotation[i]==1 || szAnnotation[i]==szArray[i], 1:length(szAnnotation))
end

function annotate!(A::AnnotatedArray{T,N},name::Symbol,annotation::SingletonArray) where {T,N}
	@assert allowedsize(size(A),size(annotation))
	A.annotations[name] = annotation
	A
end
annotate!(A::AnnotatedArray{T,N},name::Symbol,annotation::AbstractArray) where {T,N} = annotate!(A,name,SingletonArray(annotation))







size(A::AnnotatedArray) = size(A.array)
Base.IndexStyle(::Type{AnnotatedArray}) = Base.IndexStyle(Type{Array})


# retrieve as SingletonArray
annotation(A::AnnotatedArray, s::Symbol) = A.annotations[s]

# retrieve underlying annotation array
getindex(A::AnnotatedArray, s::Symbol) = A.annotations[s].array
setindex!(A::AnnotatedArray, v, s::Symbol) = annotate!(A,s,v)


getindex(A::AnnotatedArray, i::Int)	= getindex(A.array,i)
setindex!(A::AnnotatedArray, v, i::Int) = setindex(A.array,v,i)
getindex(A::AnnotatedArray, I::Int...)     = getindex(A.array,I...)
setindex!(A::AnnotatedArray, v, I::Int...) = setindex!(A.array,v,I...)


# swap annotated indices for normal indices


deannotateindex(::AnnotatedArray, dim::Int, ind) = ind # by default use the index as is

function deannotateindex(A::AnnotatedArray, dim::Int, ind::AnnotatedIndex)
	# TODO: error checking
	evaluateindex(A,ind)
end

function deannotateindex(A::AnnotatedArray, I...)
	([deannotateindex(A,i,ind) for (i,ind) in enumerate(I)]..., )
	# ([ deannotateindex(A, i, ind) for (i, ind) = enumerate(I) ]...)
	# ([ deannotateindex(A, i, ind) for (i, ind) = enumerate(I) ]...,)
end


function getindex(A::AnnotatedArray, I...)
	I = deannotateindex(A,I...)

	B = AnnotatedArray(getindex(A.array,I...))
	for (k,v) in A.annotations
		annotate!(B, k, v[I...])
	end
	B
end

function setindex!(A::AnnotatedArray, v, I...)
	I = deannotateindex(A,I...)
	setindex!(A.array,v,I...)
end



function viewimpl(A::AnnotatedArray, I...)
	I = deannotateindex(A,I...)
	B = AnnotatedArray(view(A.array,I...))
	for (k,v) in A.annotations
		annotate!(B,k,view(v,I...))
	end
	B
end


view(A::AnnotatedArray, I::Union{Int64, AbstractVector, Colon}...) = viewimpl(A,I...)
view(A::AnnotatedArray, I::Union{Int64, AbstractVector, Colon, AnnotatedIndex}...) = viewimpl(A,I...)


copy(A::AnnotatedArray) = AnnotatedArray(copy(A.array),copy(A.annotations))



# TODO: implement cat that works in multiple dimensions
#function cat(dim::Int, A::AnnotatedArray...)
function cat(dim::Int, A::AnnotatedArray, A2::AnnotatedArray...)
	A = (A,A2...) 

	# check that the arrays have compatible sizes
	nd = maximum(map(ndims,A))
	for i=1:nd
		i==dim && continue
		sz = [map( x->size(x,i), A )...] # all sizes along this dimension

		if any(sz.!=sz[1]) 
			firstdifferent = sz[sz.!=sz[1]][1]
			throw("Dimension mismatch in dimension $i (expected $(sz[1]) got $(firstdifferent))")
		end
	end

	B = AnnotatedArray(cat(dim,map(x->convert(Array,x),A)...))

	
	# handle annotations
	
	# gather all annotations that exist in at least one array
	allNames = union(map(annotationnames,A)...)

	for name in allNames
		annots = map(x->annotation(x,name), A) # get this annotation from all arrays

		
		if all( x->isequal(x,annots[1]), annots )
			# all annotations are equal

			if size(annots[1],dim) == 1
				# just keep annotation
				annotate!(B, name, copy(annots[1])) # copy is expected when concatenating
				continue
			end
		end

		# some differences between the annotations or we are concatenating along non-singleton dimension
		# TODO: handle errors better
		a = cat(dim, annots...)
		annotate!(B,name,a)
	end
	
	B
end


vcat(A::AnnotatedArray...) = cat(1,A...)
hcat(A::AnnotatedArray...) = cat(2,A...)



# # true if all annotations that are in both A and B are identical
# function compatibleannotations(A::AnnotatedArray, B::AnnotatedArray)
# 	all(n->isequal(A[n],B[n]), intersect(annotationnames(A), annotationnames(B)))
# end


# # Simple arithmetic
# function +(A::AnnotatedArray, B::AnnotatedArray)
# 	@assert size(A)==size(B)
# 	@assert compatibleannotations(A,B)

# 	R = AnnotatedArray(A.array+B.array)
# 	for n in annotationnames(A)
# 		annotate!(R, n, copy(A[n]))
# 	end
# 	for n in setdiff(annotationnames(B), annotationnames(A))
# 		annotate!(R, n, copy(B[n]))
# 	end
# 	R
# end






# # show:
# # Simple example, with one column and one row annotation
# #         CName c1  c2  c3
# # RowName       __________
# #      r1      |
# #      r2      |
# #      r3      |

# #
# # Print as 2d-slices (similar to Array printing).
# # Annotations are printed "when needed".
# # sz[1:2] = (N,1)
# # 	printed as row annotation
# # sz[1:2] = (1,M)
# # 	printed as columnt annotation
# # sz[1:2] = (N,M)
# #	printed before slice as matrix
# # sz[1:2] = (1,1)
# #	printed before the slice


# helper function for display
# TODO: there's gotta be a better way to do this
function gatherannotations(A::AnnotatedArray)
	rAnnot = Symbol[]
	cAnnot = Symbol[]
	bAnnot = Symbol[]
	sAnnot = Symbol[]

	for (k,v) in A.annotations
		if size(v,1)==1 && size(v,2)==1
			push!(sAnnot,k) # singleton annotation
		elseif size(v,1)==size(A,1) && size(v,2)==1
			push!(rAnnot,k) # row annotation
		elseif size(v,1)==1 && size(v,2)==size(A,2)
			push!(cAnnot,k) # column annotation
		elseif size(v,1)==size(A,1) && size(v,2)==size(A,2)
			push!(bAnnot,k) # both row and column annotation (i.e. matrix)
		else
			throw(InvalidStateException("Size of annotation doesn't match size of AnnotatedArray."))
		end
	end
	rAnnot, cAnnot, bAnnot, sAnnot
end

function printsingletonannotations(io::IO,A::AnnotatedArray, annot)
	for k in annot
		print(io,k,"=",A[k][1], '\t')
	end
	isempty(annot) || println(io)
end


# TODO: rewrite using smarter widths (i.e. gather width of each column first, and print using this info later)
#function showslice(io::IO,A::AnnotatedArray, rAnnot, cAnnot, bAnnot, sAnnot)
function showslice(io::IO,A::AnnotatedArray, rAnnot, cAnnot, bAnnot)
	# full annotations
	for k in bAnnot
		v = A[k]
		pre = "$k=["
		post = ']'
		indent = " "^length(pre)
		print(io,pre)
		for i=1:size(v,1)
			i==1 || print(io,indent)
			for j=1:size(v,2)
				j==1 || print(io,'\t')
				print(io,v[i,j])
			end
			i!=size(v,1) && println()
		end
		println(io,post)
	end


	# column annotations
	for k in cAnnot
		v = A[k]
		print(io,"\t"^length(rAnnot))
		print(io,k)

		for j=1:size(A,2)
			print(io,'\t', v[j])
		end
		println(io)
	end
	
	# row annotation headers
	for k in rAnnot
		v = A[k]
		print(io,k,'\t')
	end
	println(io)

	for i=1:size(A,1)
		for k in rAnnot
			v = A[k]
			print(io, v[i], '\t')
		end
		print(io,'|')
		for j=1:size(A,2)
			print(io, '\t', A[i,j])
		end
		println(io)
	end
end


function printslicedesc(io::IO,ci::CartesianIndex)
	print(io,"[:,:")
	for i in ci.I
		print(io, ',', i)
	end
	print(io,']')
end
printlnslicedesc(io::IO,ci::CartesianIndex) = (printslicedesc(io,ci); println(io))

function showall(io::IO, A::AnnotatedArray)
	# show slices
	# TODO: fix ugly code for handling annotations

	sz = size(A)

	# gather row, columnn, both, singleton annotations (from the slice perspective)
	rAnnot, cAnnot, bAnnot, sAnnot = gatherannotations(A)

	if length(sz)<=2
		printsingletonannotations(io, A, sAnnot)
		showslice(io,A,rAnnot,cAnnot,bAnnot)
	else
		trailingSz = sz[3:end]
		#nonSingletonDims = map( x->[size(annotation(A,x),(1:length(sz))...)...].>1, sAnnot )
		nonSingletonDims = Array{BitArray{1},1}(length(sAnnot))
		map!( x->[size(annotation(A,x),(1:length(sz))...)...].>1, nonSingletonDims, sAnnot )

		# show those that are singleton in all dimensions
		length(sAnnot)>0 && printsingletonannotations(io, A, sAnnot[.!map(any,nonSingletonDims)])

		prevslice = CartesianIndex((zeros(Int,length(trailingSz))...,))
		for slice in CartesianRange(trailingSz) # for each slice
			printlnslicedesc(io,slice)
			B = view(A,:,:,slice.I...) # TODO: can I avoid using CartesianIndex.I ???

			# print sAnnot that are nonsingleton in a dimension that affects the slice
			dimChanged = Bool[slice[i] != prevslice[i] for i=1:length(slice)]
			prevslice = slice

			showAnnot = Bool[any(dimChanged .& nonSingletonDims[i][3:end]) for i=1:length(sAnnot)]
			printsingletonannotations(io, B, sAnnot[showAnnot])

			showslice(io,B,rAnnot,cAnnot,bAnnot)
		end
	end


end

# simple implementation to make sure we limit the output
function show(io::IO, A::AnnotatedArray)
	sz = size(A)
	szLimits = 3*ones(Int,length(sz))
	szLimits[1] = 20
	length(szLimits)>1 && (szLimits[2] = 8)

	szSub = (min.([sz...],szLimits)...,)

	if sz==szSub
		println(io, join(sz,'x'), " ", typeof(A), ":")
		showall(A)
	else
		B = view(A, map(n->1:n, szSub)...)
		println(io, "Showing ", join(szSub,'x'), " part of ", join(sz,'x'), " ", typeof(A), ":")
		showall(B)
	end
end



# TODO: is this the proper way?
display(A::AnnotatedArray) = show(STDOUT,A)





struct Annot
	s::Symbol
end

annot(s::Symbol) = Annot(s)
annot(s::AbstractString) = Annot(symbol(s))

struct Unary{Op,T} <: AnnotatedIndex
	rhs::T
end

struct Cmp{S,Op,T} <: AnnotatedIndex
	lhs::S
	rhs::T
end



evaluatevalue(A::AnnotatedArray, a::Annot) = A[a.s][:]

evaluatevalue(A::AnnotatedArray, v::Any) = v

evaluateindex(A::AnnotatedArray, ind) = ind



# TODO: generate code for all operators using macros
evaluateindex(A::AnnotatedArray, ind::Cmp{S,:.<,T}) where {S,T} = evaluatevalue(A,ind.lhs) .< evaluatevalue(A,ind.rhs)
<(a::Annot,b::Annot) = Cmp{Annot,:.<,Annot}(a,b)
<(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:.<,T}(a,b)
<(a::T,b::Annot) where {T<:AbstractArray} = Cmp{T,:.<,Annot}(a,b)
<(a::Annot,b::T) where {T} = Cmp{Annot,:.<,T}(a,b)
<(a::T,b::Annot) where {T} = Cmp{T,:.<,Annot}(a,b)

evaluateindex(A::AnnotatedArray, ind::Cmp{S,:.<=,T}) where {S,T} = evaluatevalue(A,ind.lhs) .<= evaluatevalue(A,ind.rhs)
<=(a::Annot,b::Annot) = Cmp{Annot,:.<=,Annot}(a,b)
<=(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:.<=,T}(a,b)
<=(a::T,b::Annot) where {T<:AbstractArray} = Cmp{T,:.<=,Annot}(a,b)
<=(a::Annot,b::T) where {T} = Cmp{Annot,:.<=,T}(a,b)
<=(a::T,b::Annot) where {T} = Cmp{T,:.<=,Annot}(a,b)

evaluateindex(A::AnnotatedArray, ind::Cmp{S,:.>,T}) where {S,T} = evaluatevalue(A,ind.lhs) .> evaluatevalue(A,ind.rhs)
>(a::Annot,b::Annot) = Cmp{Annot,:.>,Annot}(a,b)
>(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:.>,T}(a,b)
>(a::T,b::Annot) where {T<:AbstractArray} = Cmp{T,:.>,Annot}(a,b)
>(a::Annot,b::T) where {T} = Cmp{Annot,:.>,T}(a,b)
>(a::T,b::Annot) where {T} = Cmp{T,:.>,Annot}(a,b)

evaluateindex(A::AnnotatedArray, ind::Cmp{S,:.>=,T}) where {S,T} = evaluatevalue(A,ind.lhs) .>= evaluatevalue(A,ind.rhs)
>=(a::Annot,b::Annot) = Cmp{Annot,:.>=,Annot}(a,b)
>=(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:.>=,T}(a,b)
>=(a::T,b::Annot) where {T<:AbstractArray} = Cmp{T,:.>=,Annot}(a,b)
>=(a::Annot,b::T) where {T} = Cmp{Annot,:.>=,T}(a,b)
>=(a::T,b::Annot) where {T} = Cmp{T,:.>=,Annot}(a,b)

evaluateindex(A::AnnotatedArray, ind::Cmp{S,:.==,T}) where {S,T} = evaluatevalue(A,ind.lhs) .== evaluatevalue(A,ind.rhs)
==(a::Annot,b::Annot) = Cmp{Annot,:.==,Annot}(a,b)
==(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:.==,T}(a,b)
==(a::T,b::Annot) where {T<:AbstractArray} = Cmp{T,:.==,Annot}(a,b)
==(a::Annot,b::T) where {T} = Cmp{Annot,:.==,T}(a,b)
==(a::T,b::Annot) where {T} = Cmp{T,:.==,Annot}(a,b)

evaluateindex(A::AnnotatedArray, ind::Cmp{S,:.!=,T}) where {S,T} = evaluatevalue(A,ind.lhs) .!= evaluatevalue(A,ind.rhs)
!=(a::Annot,b::Annot) = Cmp{Annot,:.!=,Annot}(a,b)
!=(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:.!=,T}(a,b)
!=(a::T,b::Annot) where {T<:AbstractArray} = Cmp{T,:.!=,Annot}(a,b)
!=(a::Annot,b::T) where {T} = Cmp{Annot,:.!=,T}(a,b)
!=(a::T,b::Annot) where {T} = Cmp{T,:.!=,Annot}(a,b)



# Set functionality
function evaluateindex(A::AnnotatedArray, ind::Cmp{S,:inset,T}) where {S,T}
	rhs = evaluatevalue(A,ind.rhs)
	map(x->x∈rhs, evaluatevalue(A,ind.lhs))
end
inset(a::Annot,b::T) where {T<:AbstractArray} = Cmp{Annot,:inset,T}(a,b)
inset(s::Symbol,b::T) where {T<:AbstractArray} = inset(annot(s),b)



evaluateindex(A::AnnotatedArray, ind::Cmp{S,:&,T}) where {S,T} = evaluateindex(A,ind.lhs) .& evaluateindex(A,ind.rhs)
Base. &(a::C,b::D) where {C<:Cmp,D<:Cmp} = Cmp{C,:&,D}(a,b)
Base. &(a::C,b::T) where {C<:Cmp,T}      = Cmp{C,:&,T}(a,b)
Base. &(a::T,b::C) where {C<:Cmp,T}      = Cmp{T,:&,C}(a,b)

evaluateindex(A::AnnotatedArray, ind::Cmp{S,:|,T}) where {S,T} = evaluateindex(A,ind.lhs) .| evaluateindex(A,ind.rhs)
Base. |(a::C,b::D) where {C<:Cmp,D<:Cmp} = Cmp{C,:|,D}(a,b)
Base. |(a::C,b::T) where {C<:Cmp,T}      = Cmp{C,:|,T}(a,b)
Base. |(a::T,b::C) where {C<:Cmp,T}      = Cmp{T,:|,C}(a,b)


evaluateindex(A::AnnotatedArray, ind::Unary{:!,T}) where {T}= !evaluateindex(A,ind.rhs)
!(a::C) where {C<:Cmp} = Unary{:!,C}(a)

evaluateindex(A::AnnotatedArray, ind::Unary{:~,T}) where {T}= ~evaluateindex(A,ind.rhs)
~(a::C) where {C<:Cmp} = Unary{:~,C}(a)




# function transpose(A::AnnotatedArray{T,2,S}) where {T,S}
# 	annotations = Dict{Symbol,SingletonArray}(map(x->(x[1],transpose(x[2])), A.annotations))
# 	a = AnnotatedArray(transpose(A.array))
# 	a.annotations = annotations
# 	a
# end

# function ctranspose(A::AnnotatedArray{T,2,S}) where {T,S}
# 	annotations = Dict{Symbol,SingletonArray}(map(x->(x[1],ctranspose(x[2])), A.annotations))
# 	a = AnnotatedArray(ctranspose(A.array))
# 	a.annotations = annotations
# 	a
# end



function convert(::Type{Dict}, a::AnnotatedArray)
	@assert !haskey(a.annotations,:data)
	d = Dict{Symbol,Any}()
	for (k,v) in a.annotations
		d[k] = convert(Array,v)
	end
	d[:data] = convert(Array,a)
	d
end


function convert(::Type{AnnotatedArray}, d::Dict{Symbol,Any})
	@assert haskey(d,:data)
	a = AnnotatedArray(pop!(d,:data))
	for (k,v) in d
		annotate!(a,k,v)
	end
	a
end
function convert(::Type{AnnotatedArray}, d::Dict{T,Any}) where {T<:AbstractString}
	convert(AnnotatedArray,Dict{Symbol,Any}([Symbol(k)=>v for (k,v) in d])) # can be simplified by call to map in julia 0.5	
end



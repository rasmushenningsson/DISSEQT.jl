
# internal helper stuff for handling pattern matching

SimplePattern = Union{Regex,AbstractString,Function}
Pattern = Union{SimplePattern, AbstractArray}


matchany(p::Regex, s::AbstractString) = match(p,s)!=nothing
matchany(p::AbstractString, s::AbstractString) = match(Regex(p),s)!=nothing
matchany(p::Function, s::AbstractString) = p(s)

matchany(p::AbstractArray, s::AbstractString) = any(x->match(x,s)!=nothing,p)

matchany(p::Pattern, s::AbstractArray) = map!(x->matchany(p,x),BitArray(size(s)),s)



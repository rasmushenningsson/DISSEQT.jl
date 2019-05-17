
function getreferencegenomes(syn, referenceFolder::AbstractString, referenceNames::AbstractArray)
    ids = map( x->getchildbyname(syn, referenceFolder, string(x,".fasta")), referenceNames )
    localPaths = map(x->localpath(syn,x), ids)
    referenceGenomes = map(loadfasta, localPaths) 
    
    if all(x->length(x)==1, referenceGenomes)
        referenceGenomes = map(x->x[1][2], referenceGenomes) # NB: Only for genomes with 1 segment, get the sequence part
    end

    ids, referenceGenomes
end

# flagsSet:   flags that must be set
# flagsUnset: flags that must be unset
function processbam!(bamFile::BamFile, ra; mappingQualityThreshold::Int=30, 
                     flagsSet::Int=0, flagsUnset::Int=0, ignoreChimeric::Bool=true)
    for r in bamFile
        f = flag(r)
        f&flagsSet   != flagsSet && continue
        f&flagsUnset != 0        && continue

        mapq(r)<mappingQualityThreshold && continue

        ignoreChimeric && hastag(r,tag"SA") && continue

        processread!(ra,r)
    end
end



phred2prob(q::Integer) = min(10.0^(-Float64(q)/10.0),0.75) # cap at 0.75, since more gives lower chances of this base than the others!


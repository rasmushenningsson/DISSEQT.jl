
function saveplot(pl::py.SyncPlot,format::Symbol,outPath::AbstractString; js=:remote)
    filename = string(outPath,'.',format)
    if format==:html
        py.savehtml(pl, filename, js)
    elseif format in [:png, :pdf, :svg]
        py.savefig(pl, filename)
    else
        error("Unknown image format $format")
    end
    filename
end

saveplot(pl::py.SyncPlot,formats::AbstractArray{Symbol},outPath::AbstractString;kwargs...) = map(fmt->saveplot(pl,fmt,outPath;kwargs...), formats)

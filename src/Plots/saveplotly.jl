
_savehtml(pl::py.SyncPlot,filename,js) = py.savefig(pl,filename,js=js)
function _saverendered(pl::py.SyncPlot,filename,delay)
    display(pl)
    sleep(delay) # the problem is we don't know how long
    py.savefig(pl, filename)
    close(pl)
end

function saveplot(pl::py.SyncPlot,format::Symbol,outPath::AbstractString; delay=5, js=:remote)
    filename = string(outPath,'.',format)
    if format==:html
        _savehtml(pl,filename,js)
    elseif format in [:png, :pdf, :svg]
        _saverendered(pl,filename,delay)
    else
        error("Unknown image format $format")
    end
    filename
end

saveplot(pl::py.SyncPlot,formats::AbstractArray{Symbol},outPath::AbstractString;kwargs...) = map(fmt->saveplot(pl,fmt,outPath;kwargs...), formats)



_hasbg(pl::Plot) = pl.theme.background_color != nothing
_addbg!(pl::Plot) = pl.theme.background_color=colorant"white"
_removebg!(pl::Plot) = pl.theme.background_color=nothing


# works for both plots and subplots arranged by hstack/vstack
function _saveplot(pl,format::Symbol,outPath::AbstractString; width=29.7cm, height=21cm, kwargs...)
    filename = string(outPath,'.', format==:jssvg ? "js.svg" : format)
    if format==:png
        draw(PNG(filename, width, height; kwargs...), pl)
    elseif format==:pdf
        draw(PDF(filename, width, height; kwargs...), pl)
    elseif format==:svg
        draw(SVG(filename, width, height; kwargs...), pl)
    elseif format==:jssvg
        draw(SVGJS(filename, width, height; kwargs...), pl)
    else
        error("Unknown image format $format")
    end
    filename
end


function saveplot(pl::Plot,format::Symbol,outPath::AbstractString; solidBG=format==:png, kwargs...)
    changeBG = solidBG && !_hasbg(pl)
    changeBG && _addbg!(pl)
    filename = _saveplot(pl,format,outPath; kwargs...)
    changeBG && _removebg!(pl)
    filename
end

function saveplot(pl::AbstractArray{Plot},format::Symbol,outPath::AbstractString; solidBG=format==:png, kwargs...)
    changeBG = .!map(_hasbg, pl) .& solidBG
    for (p,c) in zip(pl,changeBG)
        c && _addbg!(p)
    end

    combined = hstack(mapslices(vstack,pl,2)...) # stack vertically and then horizontally
    filename = _saveplot(combined,format,outPath; kwargs...)

    for (p,c) in zip(pl,changeBG)
        c && _removebg!(p)
    end

    filename
end

saveplot(pl,formats::AbstractArray{Symbol},outPath::AbstractString; kwargs...) = map(fmt->saveplot(pl,fmt,outPath;kwargs...), formats)

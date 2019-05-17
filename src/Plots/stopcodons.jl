
function _pointlayer(df::DataFrame, x::Symbol, y::Symbol; groupBy::Symbol=Symbol())
    groupBy==Symbol() ? layer(df,x=x,y=y,Geom.point) : layer(df,x=x,y=y,color=groupBy,Geom.point)
end
function _dataframelinearregression(df::DataFrame,x::Symbol,y::Symbol,group=[])
    dataX = convert(Array, df[x])
    dataY = convert(Array, df[y])
    xVals = collect(extrema(dataX))

    A = xVals[2]-xVals[1] > 1e-6 ?
        hcat(ones(dataX),dataX)\dataY : # linear regression: y = A[1] + A[2]x
        (0.0,0.0)
    
    isempty(group) ? DataFrame(x=xVals, y=A[1]+A[2]*xVals) :
                     DataFrame(x=xVals, y=A[1]+A[2]*xVals, group=group)
end
function _linearregressionlayer(df::DataFrame, x::Symbol, y::Symbol; groupBy::Symbol=Symbol())
    if groupBy==Symbol()
        dfLinear = _dataframelinearregression(df,x,y)
        layer(dfLinear,x=x,y=y,Geom.line)
    else
        dfLinear = DataFrame()
        for group in unique(df[groupBy])
            dfMasked = df[df[groupBy].==group, :]
            dfLinear = vcat(dfLinear, _dataframelinearregression(dfMasked,x,y,group))
        end
        layer(dfLinear,x=:x,y=:y,color=:group,Geom.line)
    end
end




function linearregressionplot(df::DataFrame, x::Symbol, y::Symbol, args...; groupBy::Symbol=Symbol(), kwargs...)
    plot( _pointlayer(df,x,y,groupBy=groupBy),
          _linearregressionlayer(df,x,y,groupBy=groupBy),
          args...;
          kwargs...)
end
function plotlabel(string::String, fontSize, plotPadding)
    plot(x=[0],y=[0],label=[string], Geom.label(position=:centered), 
         Guide.xlabel(""), Guide.ylabel(""), Guide.xticks(label=false), Guide.yticks(label=false),
         Theme(point_label_font_size=fontSize, plot_padding=plotPadding, grid_color=RGBA(1,1,1,0)))
end

function pairwisescatterplot(X::AbstractMatrix, groupBelow, groupColorsBelow, groupAbove, groupColorsAbove;
                             padding = 10.0mm/size(X,2), point_size=0.9mm, highlight_width=point_size/3)
    N = size(X,2)
    P = Matrix{Plot}(N,N)
    for j=1:N
        for i=1:N
            if i==j
                P[i,j] = plotlabel(string(i), 120mm/N, padding)
            else
                group,groupColors = i>j ? (groupBelow,groupColorsBelow) : (groupAbove,groupColorsAbove)
                x,y = X[:,i],-X[:,j]
                df = DataFrame(x=x, y=y, group=group)
                P[i,j] = plot(df, x=:x, y=:y, color=:group, groupColors, 
                              Guide.xlabel(""), Guide.ylabel(""),
                              Guide.xticks(ticks=[extrema(x)...],label=false), Guide.yticks(ticks=[extrema(y)...],label=false),
                              Theme(key_position=:none,  # disable legend
                                    plot_padding=padding,
                                    point_size=point_size,
                                    highlight_width=highlight_width))
            end
        end
    end
    P
end
pairwisescatterplot(X, group, groupColors) = pairwisescatterplot(X, group, groupColors, group, groupColors)

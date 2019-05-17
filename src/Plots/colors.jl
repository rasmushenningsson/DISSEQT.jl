

# creates a Scale.color_discrete_manual based on a vector / DataFrame column and a Dict
# every unique value in the vector must be a key in the dict
# returns Scale.color_discrete_manual with color order such that Gadfly plotting with color=myColumn will give colors as specified by the Dict.
groupcolors(v, dict::Dict) = Scale.color_discrete_manual(map(x->dict[x], unique(va))...)


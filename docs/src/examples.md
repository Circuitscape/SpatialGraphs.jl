# Examples

## Habitat Corridor Mapping

In this example, we'll map a habitat corridor in central Maryland. We'll use
SpatialGraphs.jl to convert the landscape into a graph, then calculate cost 
distances, and finally convert the result back into a raster to display a
least cost corridor in space.

First, install the necessary packages and import them:

```julia
using Pkg
Pkg.add(["Rasters", "Plots", "Graphs", "SimpleWeightedGraphs", "SpatialGraphs"])
using SpatialGraphs, Rasters, Plots, Graphs, SimpleWeightedGraphs
```

```@setup corridors
using Pkg
Pkg.add(["Rasters", "Graphs", "SimpleWeightedGraphs"])
using SpatialGraphs, Rasters, Graphs, SimpleWeightedGraphs

url_base = "https://raw.githubusercontent.com/Circuitscape/datasets/main/"
download(string(url_base, "data/nlcd_2016_frederick_md.tif"),
         "nlcd_2016_frederick_md.tif")

land_cover = Raster("nlcd_2016_frederick_md.tif")
nothing
```

Now, let's download the land cover map that we'll use as the basis for defining
the weights of our graph, and plot it. Graph weights correspond to the cost of 
an animal moving from one graph vertex to another. This is also commonly 
referred to as "resistance" in the context of connectivity modeling.

```julia
url_base = "https://raw.githubusercontent.com/Circuitscape/datasets/main/"
download(string(url_base, "nlcd_2016_frederick_md.tif"),
         "nlcd_2016_frederick_md.tif")
values = [11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95]
palette = ["#476BA0", "#DDC9C9", "#D89382", "#ED0000", "#AA0000",
           "#b2b2b2", "#68AA63", "#1C6330", "#B5C98E", "#CCBA7C",
           "#E2E2C1", "#DBD83D", "#AA7028", "#BAD8EA", "#70A3BA"]

plot(Raster("nlcd_2016_frederick_md.tif"), title = "Land Cover Type", 
     xlabel = "Easting", ylabel = "Northing", 
     seriescolor = cgrad(palette, (values .- 12) ./ 84,
     categorical = true), size = (700, 640))
```
```@raw html
<img src='../figs/mdlc.png' width=500><br>
```

Load the land cover data as a Raster, read it into memory, then convert it to 
Float64.

```@example corridors
landcover = convert.(
    Float64,
    read(
        Raster(
            "nlcd_2016_frederick_md.tif",
            missingval = -9999.0
        )
    )
)

```

Now, we need to convert land cover into resistance by reclassifying the values
in the land cover raster to their corresponding costs/weights.

```@example corridors
# Copy the landcover object to use as a place holder for resistance
resistance = deepcopy(landcover)

# Define how values should be reclassified. Original values are on the left,
# and new values are on the right.
reclass_table = [
    11.	100; # Water
    21	500; # Developed, open space
    22	1000; # Developed, low intensity
    23	-9999; # Developed, medium intensity
    24	-9999; # Developed, high intensity
    31	100; # Barren land
    41	1; # Deciduous forest
    42	1; # Evergreen forest
    43	1; # Mixed forest
    52	20; # Shrub/scrub
    71	30; # Grassland/herbaceous
    81	200; # Pasture/hay
    82	300; # Cultivated crops
    90	20; # Woody wetlands
    95	30; # Emergent herbaceous wetlands
]

# Redefine the values in resistance based on the corresponding land cover
# values, using reclass_table 
for i in 1:(size(reclass_table)[1])
    resistance[landcover .== reclass_table[i, 1]] .= reclass_table[i, 2]
end
```

Now that we have a resistance raster, we can start using SpatialGraphs.jl!

Convert the resistance raster into a [`WeightedRasterGraph`](@ref WeightedRasterGraph):

```@example corridors
wrg = weightedrastergraph(resistance)
```

Now that we have a `WeightedRasterGraph`, we can make use of the functions in
Graphs.jl.

Let's calculate cost distances from two different vertices. Vertex IDs, 
the second argument to `dijkstra_shortest_paths`, correspond to the values in 
wrg.vertex_raster.

First, we'll use `dijkstra_shortest_paths` to calculate to total cost distance
from each vertex of interest to all other vertices in the graph.

```@example corridors
cd_1 = dijkstra_shortest_paths(wrg, 1)
cd_2 = dijkstra_shortest_paths(wrg, 348021)
```

Now, create placeholder rasters to store the cost distance values for each 
vertex, and populate them using the cost distances calculated
above. Finally, we'll sum the rasters to enable us to delineate a corridor. 
The value at each pixel corresponds to the total cost distance of the least cost
path between vertices 1 and 348021 that intersects that pixel.

```@example corridors
# Create rasters
cd_1_raster = Raster( # Start with an empty raster that will be populated
    fill(float(resistance.missingval), size(wrg.vertex_raster)),
    dims(resistance),
    missingval=resistance.missingval
)
cd_2_raster = Raster(
    fill(float(wrg.vertex_raster.missingval), size(wrg.vertex_raster)),
    dims(resistance),
    missingval=resistance.missingval
)

# Populate with cost distances
cd_1_raster[wrg.vertex_raster .!= wrg.vertex_raster.missingval] .= cd_1.dists[
    wrg.vertex_raster[wrg.vertex_raster.!= wrg.vertex_raster.missingval]
]
cd_2_raster[wrg.vertex_raster .!= wrg.vertex_raster.missingval] .= cd_2.dists[
    wrg.vertex_raster[wrg.vertex_raster.!= wrg.vertex_raster.missingval]
]

# Sum the rasters to enalbe corridor mapping
cwd_sum = cd_1_raster + cd_2_raster
```

Setting some land cover classes to NoData resulted in some vertices being
completely isolated from the vertices from which we calculated cost distances
above. In these cases, `dijkstra_shortest_paths` returns a value of `Inf` for
the cost distance. Let's get rid of those by setting them to the NoData value
of our raster.

```@example corridors
cwd_sum[cwd_sum .== Inf] .= cwd_sum.missingval
```

Finally, let's plot the result! We'll use a threshold to set pixels above
a certain value to NoData, so we will only retain the pixels within our
"corridor".

```julia
corridor = copy(cwd_sum)
corridor[cwd_sum .>= (minimum(cwd_sum[cwd_sum .>= 0]) + 250)] .= -9999

plot(
    cwd_sum,
    title = "Cost Distance Corridor",
    xlabel = "Easting", ylabel = "Northing",
    colorbar_title = " \nCost Distance Sum",
    size = (650, 475)
)
```

```@raw html
<img src='../figs/corridor_plot.png' width=500><br>
```

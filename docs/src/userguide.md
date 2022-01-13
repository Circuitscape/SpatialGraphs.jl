# User Guide

## Building Graphs from Rasters

SpatialGraphs.jl offers several functions for constructing graphs from raster
data, which are detailed below. Once you have converted your data to an
[`AbstractSpatialGraph`](@ref AbstractSpatialGraph), you can analyze the graph
using funtions from Graphs.jl. Following you analysis, you can then use the 
spatial information stored in the `AbstractSpatialGraph` to map values 
(such as betweenness or cost distance, for example) back to the appropriate 
points in space. See the [Examples](@ref Examples) section for a detailed 
demonstration of how this can be done.

### Simple Graphs
```@docs
rastergraph
make_simple_raster_graph
```

### Weighted Graphs
```@docs
weightedrastergraph
make_weighted_raster_graph
```

### Making vertex rasters
In SpatialGraphs.jl, vertex rasters serve as the spatial reference for the 
vertices in a graph. SpatialGraphs.jl provides functions to generate these
rasters. Often, it is done internally and the end user doesn't need to use these
functions, but there are some cases where it will be necessary. For example, it 
is often the case that a user may want a single graph vertex to occupy multiple
pixels in space (e.g. when modeling habitat connectivity between two
protected areas). SpatialGraphs.jl enables this by offering a method for the 
`make_vertex_raster` function (below) that accepts a raster
representing these patches as an argument.
```@docs
make_vertex_raster
```
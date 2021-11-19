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
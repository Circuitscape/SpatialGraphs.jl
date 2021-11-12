# [Graph Types in SpatialGraphs.jl](@id graph_types)

At this time, only raster-based graph types have been developed (and 
vertex locations are stored in a `Rasters.Raster`), but there are plans to 
eventually implement graph types for vector data as well.

## Abstract Types
The overarching type in SpatialGraphs.jl is the `AbstractSpatialGraph`. All
`AbstractSpatialGraph` subtypes contain a field called `graph`, which is itself
an `AbstractGraph`.
```@docs
AbstractSpatialGraph
```

The `AbstractRasterGraph` type is a subtype of `AbstractSpatialGraph`. All
`AbstractRasterGraph` subtypes contain a field called `vertex_raster`, which is 
a `Rasters.Raster` that describes the spatial locations for the vertices in 
the graph. If your graph has 10 vertices, then the corresponding `vertex_raster` 
must have 10 unique values, starting at 1 (which corresponds to the first vertex
in the graph) and going up to 10 (which corresponds to the tenth vertex in the 
graph). To find the location of a given graph vertex, _n_, you simply identify 
the pixel(s) with a value of _n_. For pixels/elements in `vertex_raster` for
which there shouldn't be a graph vertex, use a value of 0.
```@docs
AbstractRasterGraph
```

## AbstractRasterGraph Subtypes
SpatialGraphs.jl has raster graph types for undirected, directed, unweighted,
and weighted graphs, each detailed below.

```@docs
RasterGraph
RasterDiGraph
WeightedRasterGraph
WeightedRasterDiGraph
```

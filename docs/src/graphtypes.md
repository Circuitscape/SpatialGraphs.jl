# Graph Types in SpatialGraphs.jl

At this time, only raster-based graph types have been developed (and 
vertex locations are stored in a `GeoData.GeoArray`), but there are plans to 
eventually implement graph types for vector data as well.

## Abstract Types
The overarching type in SpatialGraphs.jl is the `AbstractSpatialGraph`. All
`AbstractSpatialGraph` subtypes contain a field called `graph`, which is itself
an `AbstractGraph`.
```@docs
AbstractSpatialGraph
```

The `AbstractRasterGraph` type is as subtype of `AbstractSpatialGraph`. All
`AbstractRasterGraph` subtypes contain a field called `vertex_raster`, which is 
a `GeoData.GeoArray` that describes the spatial locations for the vertices in 
the graph. If your graph has 10 vertices, then the corresponding `GeoArray` 
must have 10 unique values, starting at 1 (which corresponds to the first vertex
in the graph) and going up to 10 (which corresponds to the tenth vertex in the 
graph). To find the location of a given graph vertex, _n_, you simply identify 
the pixel(s) with a value of _n_.
```@docs
AbstractRasterGraph
```

## AbstractRasterGraph Subtypes
SpatialGraphs.jl has raster graph types for undirected, directed, unweighted,
and weighted graphsm each detailed below.

```@docs
SimpleRasterGraph
SimpleRasterDiGraph
WeightedRasterGraph
WeightedRasterDiGraph
```

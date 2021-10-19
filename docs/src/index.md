# SpatialGraphs.jl

SpatialGraphs.jl introduces the `AbstractSpatialGraph`. `AbstractSpatialGraphs` 
are a subtype of `Graphs.AbstractGraph`, and can be weighted or directed.
SpatialGraphs.jl is useful for turning spatial data into graphs. This can be
useful for landscape connectivity analysis, hydrology, and other spatial
network processes. AbstractSpatialGraphs are AbstractGraphs, so methods from 
Graphs.jl work right out of the box. Go to [Graph Types](@ref graph_types) 
for more details on the graph types implemented in this package.


# SpatialGraphs.jl

SpatialGraphs.jl is  low level package that offers types and constructors for
spatially referenced graphs. The primary type is `AbstractSpatialGraph`, which is a 
`LightGraphs.AbstractGraph`. `AbstractSpatialGraphs`s work out of the box with 
functions from LightGraphs.jl and, where applicable, SimpleWeightedGraphs.jl.
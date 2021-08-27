# SpatialGraphs.jl

[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.circuitscape.org/SpatialGraphs.jl/latest)

SpatialGraphs.jl introduces the `AbstractSpatialGraph`. `AbstractSpatialGraphs` 
are a subtype of `LightGraphs.AbstractGraph`, and can be weighted or directed. 
AbstractSpatialGraphs are AbstractGraphs, methods from LightGraphs.jl work right
out of the box.

`AbstractSpatialGraph`s themselves contain an `AbstractGraph` in addition to 
metadata that details the spatial location of each vertex in the
graph. At this time, only raster-based graph types have been developed (and 
vertex locations are stored in a `GeoData.GeoArray`), but there are plans to 
implement graph types for vector data as well.

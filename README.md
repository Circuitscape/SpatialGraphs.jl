# SpatialGraphs.jl

_Types and constructors for spatially referenced graphs._

| **Documentation** | **Build Status**|
|:-----------------------------------------------------:|:------------------------------------:|
[![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.circuitscape.org/SpatialGraphs.jl/latest) | [![Build Status](https://github.com/Circuitscape/SpatialGraphs.jl/workflows/CI/badge.svg)](https://github.com/Circuitscape/SpatialGraphs.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/Circuitscape/SpatialGraphs.jl/branch/main/graph/badge.svg?token=67OX4UPWOL)](https://codecov.io/gh/Circuitscape/SpatialGraphs.jl)

SpatialGraphs.jl introduces the `AbstractSpatialGraph` and other subtypes. The `AbstractSpatialGraph` 
is a `Graphs.AbstractGraph`, and can be weighted or directed. SpatialGraphs.jl extends the
Graphs interface, so methods from Graphs.jl work right out of the box.

An `AbstractSpatialGraph` is a composite type that contains an `AbstractGraph` in addition to 
metadata that details the spatial location of each vertex in the
graph. At this time, only raster-based graph types have been developed (and 
vertex locations are stored in a `Rasters.Raster`), but there are plans to 
implement graph types for vector data as well.

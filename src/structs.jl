## All types and structs for SpatialGraphs.jl
"""
    AbstractSpatialGraph{T}

An abstract type representing a spatially referenced graph.

An AbstractSpatialGraph must contain the following attributes/element:

- `graph::AbstractGraph`
"""
abstract type AbstractSpatialGraph{T} <: AbstractGraph{T} end

"""
    AbstractRsterGraph{T}

An abstract type representing a spatially referenced graph, with graph vertices corresponding to pixels in a raster.

An AbstractSpatialGraph must contain the following attributes/element:

- `graph::AbstractGraph`
- `node_raster::GeoArray`
"""
abstract type AbstractRasterGraph{T} <: AbstractSpatialGraph{T} end

"""
    WeightedRasterGraph{T}

A composite type for a spatially referenced weighted graph. Vertices are spatially referenced based on a raster.
"""
mutable struct WeightedRasterGraph{T<:Integer, U<:Real} <: AbstractRasterGraph{T}
    graph::SimpleWeightedGraph{T, U} # a SimpleWeightedGraph with edge weights
    node_raster::GeoArray # A GeoArray raster, where a pixel's value denotes its vertex ID in the graph
end
"""
    WeightedRasterDiGraph{T}

A composite type for a spatially referenced, weighted, directed graph. Vertices are spatially referenced based on a raster.
"""
mutable struct WeightedRasterDiGraph{T<:Integer, U<:Real} <: AbstractRasterGraph{T}
    graph::SimpleWeightedDiGraph{T, U} # a SimpleWeightedDiGraph with edge weights
    node_raster::GeoArray # A GeoArray raster, where a pixel's value denotes its vertex ID in the graph
end

"""
    SimpleRasterGraph{T}

A composite type for a spatially referenced graph. Vertices are spatially referenced based on a raster.

A RasterGraph
"""
mutable struct SimpleRasterGraph{T<:Integer} <: AbstractRasterGraph{T}
    graph::SimpleGraph{T} # A SimpleGraph
    node_raster::GeoArray # A GeoArray raster, where a pixel's value denotes its vertex ID in the graph
end

"""
    SimpleRasterDiGraph{T}

A composite type for a spatially referenced directed graph. Vertices are spatially referenced based on a raster.
"""
mutable struct SimpleRasterDiGraph{T<:Integer} <: AbstractRasterGraph{T}
    graph::SimpleDiGraph{T} # A SimpleDiGraph
    node_raster::GeoArray # A GeoArray raster, where a pixel's value denotes its vertex ID in the graph
end

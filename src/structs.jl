## All types and structs for SpatialGraphs.jl
"""
    AbstractSpatialGraph{T}

An abstract type representing a spatially referenced graph.

An AbstractSpatialGraph must contain the following field:

- `graph::AbstractGraph`
"""
abstract type AbstractSpatialGraph{T} <: AbstractGraph{T} end

"""
    AbstractRasterGraph{T}

An abstract type representing a spatially referenced graph, with graph vertices corresponding to pixels in a raster.

An AbstractRasterGraph must contain the following fields:

- `graph::AbstractGraph`
- `vertex_raster::Raster`
"""
abstract type AbstractRasterGraph{T} <: AbstractSpatialGraph{T} end

"""
    WeightedRasterGraph{T}

A composite type for a spatially referenced weighted graph. Vertices are spatially referenced based on a raster.
"""
mutable struct WeightedRasterGraph{T<:Integer, U<:Real} <: AbstractRasterGraph{T}
    graph::SimpleWeightedGraph{T, U} # a SimpleWeightedGraph with edge weights
    vertex_raster::Raster # A Raster raster, where a pixel's value denotes its vertex ID in the graph
end
"""
    WeightedRasterDiGraph{T}

A composite type for a spatially referenced, weighted, directed graph. Vertices are spatially referenced based on a raster.
"""
mutable struct WeightedRasterDiGraph{T<:Integer, U<:Real} <: AbstractRasterGraph{T}
    graph::SimpleWeightedDiGraph{T, U} # a SimpleWeightedDiGraph with edge weights
    vertex_raster::Raster # A Raster raster, where a pixel's value denotes its vertex ID in the graph
end

"""
    RasterGraph{T}

A composite type for a spatially referenced graph. Vertices are spatially referenced based on a raster.
"""
mutable struct RasterGraph{T<:Integer} <: AbstractRasterGraph{T}
    graph::Graph{T} # A Graph
    vertex_raster::Raster # A Raster raster, where a pixel's value denotes its vertex ID in the graph
end

"""
    RasterDiGraph{T}

A composite type for a spatially referenced directed graph. Vertices are spatially referenced based on a raster.
"""
mutable struct RasterDiGraph{T<:Integer} <: AbstractRasterGraph{T}
    graph::DiGraph{T} # A DiGraph
    vertex_raster::Raster # A Raster raster, where a pixel's value denotes its vertex ID in the graph
end

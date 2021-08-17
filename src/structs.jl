# An AbstractSpatialGraph must contain the following elements:
# - graph::AbstractGraph
abstract type AbstractSpatialGraph{T} <: AbstractGraph{T} end
abstract type AbstractSpatialWeightedGraph{T} <: AbstractSpatialGraph{T} end

mutable struct WeightedRasterGraph{T<:Integer, U<:Real} <: AbstractSpatialWeightedGraph{T}
    graph::SimpleWeightedGraph{T, U} # a SimpleWeightedGraph with edge weights
    node_array::Matrix{T} # An array of nodes, where index corresponds to position in geographic space
    wkt::String # A WellKnownText representation of the original raster's projection, from ArchGDAL.getproj()
    transform::Vector # The geo transform, from ArchGDAL.getgeotransform()
end

mutable struct WeightedRasterDiGraph{T<:Integer, U<:Real} <: AbstractSpatialWeightedGraph{T}
    graph::SimpleWeightedDiGraph{T, U} # a SimpleWeightedDiGraph with edge weights
    node_array::Matrix{T} # An array of nodes, where index corresponds to position in geographic space
    wkt::String # A WellKnownText representation of the original raster's projection, from ArchGDAL.getproj()
    transform::Vector # The geo transform, from ArchGDAL.getgeotransform()
end

mutable struct RasterGraph{T<:Integer} <: AbstractSpatialGraph{T}
    graph::SimpleGraph{T} # A SimpleGraph
    node_matrix::Matrix{T} # An array of nodes, where index corresponds to position in geographic space
    wkt::String # A WellKnownText representation of the original raster's projection, from ArchGDAL.getproj()
    transform::Vector # The geo transform, from ArchGDAL.getgeotransform()
end

mutable struct RasterDiGraph{T<:Integer} <: AbstractSpatialGraph{T}
    graph::SimpleDiGraph{T} # A SimpleDiGraph
    node_array::Matrix{T} # An array of nodes, where index corresponds to position in geographic space
    wkt::String # A WellKnownText representation of the original raster's projection, from ArchGDAL.getproj()
    transform::Vector # The geo transform, from ArchGDAL.getgeotransform()
end

### LightGraphs interface
nv(g::AbstractSpatialGraph) = nv(g.graph)
ne(g::AbstractSpatialGraph) = ne(g.graph)
vertices(g::AbstractSpatialGraph) = vertices(g.graph)
edges(g::AbstractSpatialGraph) = edges(g.graph)
Base.eltype(g::AbstractSpatialGraph) = eltype(g.graph)
edgetype(g::AbstractSpatialGraph) = edgetype(g.graph)
has_edge(g::AbstractSpatialGraph, s, d) = has_edge(g.graph, s, d)
has_vertex(g::AbstractSpatialGraph, v) = has_vertex(g.graph, v)
inneighbors(g::AbstractSpatialGraph, v) = inneighbors(g.graph, v)
outneighbors(g::AbstractSpatialGraph, v) = outneighbors(g.graph, v)
is_directed(g::AbstractSpatialGraph) = is_directed(g.graph)
Base.zero(g::AbsractSpatialGraph) = zero(g.graph)

### SimpleWeightedGraphs interface
get_weight(g::AbstractWeightedSpatialGraph, u::Integer, v::Integer) = get_weight(g.graph, u, v)

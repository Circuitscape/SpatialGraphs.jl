## Methods for the Graphs ad SimpleWeightedGraphs interfaces
## With these methods defined, functions from Graphs should "just work"
import Graphs: 
    nv, ne, vertices, edges, eltype, edgetype, has_edge, has_vertex,
    inneighbors, outneighbors, is_directed, add_edge!

import SimpleWeightedGraphs:
    add_edge!, get_weight, add_vertex!, vertices

import Base: zero

### Graphs interface
nv(g::AbstractSpatialGraph) = nv(g.graph)
ne(g::AbstractSpatialGraph) = ne(g.graph)
vertices(g::AbstractSpatialGraph) = vertices(g.graph)
edges(g::AbstractSpatialGraph) = edges(g.graph)
eltype(g::AbstractSpatialGraph) = eltype(g.graph)
edgetype(g::AbstractSpatialGraph) = edgetype(g.graph)
has_edge(g::AbstractSpatialGraph, s, d) = has_edge(g.graph, s, d)
has_vertex(g::AbstractSpatialGraph, v) = has_vertex(g.graph, v)
inneighbors(g::AbstractSpatialGraph, v) = inneighbors(g.graph, v)
outneighbors(g::AbstractSpatialGraph, v) = outneighbors(g.graph, v)
is_directed(g::AbstractSpatialGraph) = is_directed(g.graph)
function Base.zero(g::AbstractRasterGraph)
    if g.graph isa Graph
        RasterGraph(
            zero(typeof(g.graph)),
            g.vertex_raster
        )
    elseif g.graph isa DiGraph
        RasterDiGraph(
            zero(typeof(g.graph)),
            g.vertex_raster
        )
    elseif g.graph isa SimpleWeightedGraph
        WeightedRasterGraph(
            SimpleWeightedGraph{eltype(g.graph), weighttype(g.graph)}(),
            g.vertex_raster
        )
    elseif g.graph isa SimpleWeightedDiGraph
        WeightedRasterDiGraph(
            SimpleWeightedDiGraph{eltype(g.graph), weighttype(g.graph)}(),
            g.vertex_raster
        )
    end
end
add_edge!(g::AbstractSpatialGraph, a::Integer, b::Integer, c::Number) = add_edge!(g.graph, a, b, c)
add_vertex!(g::AbstractSpatialGraph) = add_vertex!(g.graph)

### SimpleWeightedGraphs
get_weight(g::WeightedRasterGraph, a::Integer, b::Integer) = get_weight(g.graph, a, b)
get_weight(g::WeightedRasterDiGraph, a::Integer, b::Integer) = get_weight(g.graph, a, b)
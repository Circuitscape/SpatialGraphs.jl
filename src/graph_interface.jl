## Methods for the Graphs ad SimpleWeightedGraphs interfaces
## With these methods defined, functions from Graphs should "just work"
import Graphs: 
    nv, ne, vertices, edges, eltype, edgetype, has_edge, has_vertex,
    inneighbors, outneighbors, is_directed, add_edge!

import SimpleWeightedGraphs:
    add_edge!, get_weight, add_vertex!

import Base: zero

### Graphs interface
function nv(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.nv(g.graph)
    else
        nv(g.graph)
    end
end

function ne(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.ne(g.graph)
    else
        ne(g.graph)
    end
end

function vertices(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.vertices(g.graph)
    else
        vertices(g.graph)
    end
end

function edges(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.edges(g.graph)
    else
        edges(g.graph)
    end
end

function eltype(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.eltype(g.graph)
    else
        eltype(g.graph)
    end
end

function edgetype(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.edgetype(g.graph)
    else
        edgetype(g.graph)
    end
end

function has_edge(g::AbstractSpatialGraph, s, d)
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.has_edge(g.graph, s, d)
    else
        has_edge(g.graph, s, d)
    end    
end

function has_vertex(g::AbstractSpatialGraph, v)
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.has_vertex(g.graph, v)
    else
        has_vertex(g.graph, v)
    end    
end

function inneighbors(g::AbstractSpatialGraph, v)
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.inneighbors(g.graph, v)
    else
        inneighbors(g.graph, v)
    end    
end

function outneighbors(g::AbstractSpatialGraph, v)
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.outneighbors(g.graph, v)
    else
        outneighbors(g.graph, v)
    end    
end

function is_directed(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.is_directed(g.graph)
    else
        is_directed(g.graph)
    end
end

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

function add_edge!(g::AbstractSpatialGraph, a::Integer, b::Integer, c::Number)
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.add_edge!(g.graph, a, b, c)
    else
        add_edge!(g.graph, a, b, c)
    end
end

function add_vertex!(g::AbstractSpatialGraph) 
    if typeof(g.graph) <: Union{SimpleWeightedGraph, SimpleWeightedDiGraph}
        SimpleWeightedGraphs.add_vertex!(g.graph)
    else
        add_vertex!(g.graph)
    end
end

### SimpleWeightedGraphs
get_weight(g::WeightedRasterGraph, a::Integer, b::Integer) = get_weight(g.graph, a, b)
get_weight(g::WeightedRasterDiGraph, a::Integer, b::Integer) = get_weight(g.graph, a, b)
## Methods for the LightGraphs ad SimpleWeightedGraphs interfaces
## With these methods defined, functions from LightGraphs should "just work"

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

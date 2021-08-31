module SpatialGraphs

using LightGraphs, SimpleWeightedGraphs, SparseArrays, GeoData

include("structs.jl")
include("graph_interface.jl")
include("rastergraphs.jl")
include("utils.jl")
include("show.jl")

## Types and Structs
export AbstractSpatialGraph, AbstractRasterGraph, SimpleRasterGraph,
       SimpleRasterDiGraph, WeightedRasterGraph, WeightedRasterDiGraph

## Raster graph
export make_weighted_raster_graph, weightedrastergraph,
    make_simple_raster_graph, simplerastergraph


end # module

module SpatialGraphs

using Graphs, SimpleWeightedGraphs, SparseArrays, Rasters

include("structs.jl")
include("graph_interface.jl")
include("rastergraphs.jl")
include("utils.jl")
include("show.jl")

## Types and Structs
export AbstractSpatialGraph, AbstractRasterGraph, RasterGraph,
       RasterDiGraph, WeightedRasterGraph, WeightedRasterDiGraph

## Raster graph
export make_weighted_raster_graph, weightedrastergraph,
    make_simple_raster_graph, rastergraph

end # module

module SpatialGraphs

using LightGraphs, SimpleWeightedGraphs, ArchGDAL, GeoData

include("structs.jl")
include("graph_interface.jl")
include("rastergraphs.jl")
include("utils.jl")
## Types and Structs
export AbstractSpatialGraph, AbstractRasterGraph, SimpleRasterGraph,
       SimpleRasterDiGraph, WeightedRasterGraph, WeightedRasterDiGraph

## Raster graph
export make_raster_graph, weightedrastergraph

end # module

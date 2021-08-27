module SpatialGraphs

using LightGraphs, SimpleWeightedGraphs, ArchGDAL

include("structs.jl")
include("graph_interface.jl")
## Types and Structs
export AbstractSpatialGraph, AbstractSpatialWeightedGraph, RasterGraph,
       RasterDiGraph, WeightedRasterGraph, WeightedRasterDiGraph

end # module

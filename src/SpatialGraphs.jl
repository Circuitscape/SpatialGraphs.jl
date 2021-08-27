module SpatialGraphs

using LightGraphs, SimpleWeightedGraphs, ArchGDAL, GeoData

include("structs.jl")
include("graph_interface.jl")
## Types and Structs
export AbstractSpatialGraph, AbstractRasterGraph, SimpleRasterGraphs,
       SimpleRasterDiGraph, WeightedRasterGraph, WeightedRasterDiGraph

end # module

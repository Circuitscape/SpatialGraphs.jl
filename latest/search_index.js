var documenterSearchIndex = {"docs":
[{"location":"graphtypes/#graph_types","page":"Graph Types","title":"Graph Types in SpatialGraphs.jl","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"At this time, only raster-based graph types have been developed (and  vertex locations are stored in a Rasters.Raster), but there are plans to  eventually implement graph types for vector data as well.","category":"page"},{"location":"graphtypes/#Abstract-Types","page":"Graph Types","title":"Abstract Types","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"The overarching type in SpatialGraphs.jl is the AbstractSpatialGraph. All AbstractSpatialGraph subtypes contain a field called graph, which is itself an AbstractGraph.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"AbstractSpatialGraph","category":"page"},{"location":"graphtypes/#SpatialGraphs.AbstractSpatialGraph","page":"Graph Types","title":"SpatialGraphs.AbstractSpatialGraph","text":"AbstractSpatialGraph{T}\n\nAn abstract type representing a spatially referenced graph.\n\nAn AbstractSpatialGraph must contain the following field:\n\ngraph::AbstractGraph\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"The AbstractRasterGraph type is a subtype of AbstractSpatialGraph. All AbstractRasterGraph subtypes contain a field called vertex_raster, which is  a Rasters.Raster that describes the spatial locations for the vertices in  the graph. If your graph has 10 vertices, then the corresponding vertex_raster  must have 10 unique values, starting at 1 (which corresponds to the first vertex in the graph) and going up to 10 (which corresponds to the tenth vertex in the  graph). To find the location of a given graph vertex, n, you simply identify  the pixel(s) with a value of n. For pixels/elements in vertex_raster for which there shouldn't be a graph vertex, use a value of 0.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"AbstractRasterGraph","category":"page"},{"location":"graphtypes/#SpatialGraphs.AbstractRasterGraph","page":"Graph Types","title":"SpatialGraphs.AbstractRasterGraph","text":"AbstractRasterGraph{T}\n\nAn abstract type representing a spatially referenced graph, with graph vertices corresponding to pixels in a raster.\n\nAn AbstractRasterGraph must contain the following fields:\n\ngraph::AbstractGraph\nvertex_raster::Raster\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#AbstractRasterGraph-Subtypes","page":"Graph Types","title":"AbstractRasterGraph Subtypes","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"SpatialGraphs.jl has raster graph types for undirected, directed, unweighted, and weighted graphs, each detailed below.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"RasterGraph\nRasterDiGraph\nWeightedRasterGraph\nWeightedRasterDiGraph","category":"page"},{"location":"graphtypes/#SpatialGraphs.RasterGraph","page":"Graph Types","title":"SpatialGraphs.RasterGraph","text":"RasterGraph{T}\n\nA composite type for a spatially referenced graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#SpatialGraphs.RasterDiGraph","page":"Graph Types","title":"SpatialGraphs.RasterDiGraph","text":"RasterDiGraph{T}\n\nA composite type for a spatially referenced directed graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#SpatialGraphs.WeightedRasterGraph","page":"Graph Types","title":"SpatialGraphs.WeightedRasterGraph","text":"WeightedRasterGraph{T}\n\nA composite type for a spatially referenced weighted graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#SpatialGraphs.WeightedRasterDiGraph","page":"Graph Types","title":"SpatialGraphs.WeightedRasterDiGraph","text":"WeightedRasterDiGraph{T}\n\nA composite type for a spatially referenced, weighted, directed graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"userguide/#User-Guide","page":"User Guide","title":"User Guide","text":"","category":"section"},{"location":"userguide/#Building-Graphs-from-Rasters","page":"User Guide","title":"Building Graphs from Rasters","text":"","category":"section"},{"location":"userguide/","page":"User Guide","title":"User Guide","text":"SpatialGraphs.jl offers several functions for constructing graphs from raster data.","category":"page"},{"location":"userguide/#Simple-Graphs","page":"User Guide","title":"Simple Graphs","text":"","category":"section"},{"location":"userguide/","page":"User Guide","title":"User Guide","text":"rastergraph\nmake_simple_raster_graph","category":"page"},{"location":"userguide/#SpatialGraphs.rastergraph","page":"User Guide","title":"SpatialGraphs.rastergraph","text":"RasterGraph(\n    raster::Raster;\n    directed::Bool = true,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false\n)\n\nConstruct a RasterGraph or RasterDiGraph (if  directed = true) from a raster dataset.\n\nParameters\n\nraster: A Rasters.Raster on which to base the graph. Any pixel in raster  with a value not equal to raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). The values in the raster can also be used to determine which vertices to connect. See condition below for more information.\n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition: A function that compares the values in raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations that have a lower value in raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use > for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (raster.missingval). Note that if using an inequality function (or any function where the result depends on argument position), directed  should be set to true. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\n\n\n\n\n","category":"function"},{"location":"userguide/#SpatialGraphs.make_simple_raster_graph","page":"User Guide","title":"SpatialGraphs.make_simple_raster_graph","text":"make_simple_raster_graph(\n    raster::Raster,\n    vertex_raster::Raster;\n    directed::Bool = false,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false,\n    combine::Function = min\n)\n\nConstruct a SimpleGraph or SimpleDiGraph (if directed = true) based on  two rasters (raster and vertex_raster). This function is useful  if you already have a custom vertex raster and don't want SpatialGraphs.jl to make one for you. The vertex raster denotes the spatial locations of each vertex in the graph, and raster is used to construct the graph and determine which vertices to connect.\n\nParameters\n\nraster: A Rasters.Raster on which to base the graph. Any pixel in raster  with a value not equal to raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). The values in the raster can also be used to determine which vertices to connect. See condition below for more information.\n\nvertex_raster: A Rasters.Raster with integer values ranging from 1:n,  where n is the number of unique vertices in the graph. \n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition: A function that compares the values in condition_raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations with a lower value in condition_raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use < for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (condition_raster.missingval) in condition_raster. Note that if using an inequality function (or any function where the result depends on argument  position), the graph must be directed. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in condition_raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\nconnect_using_avg_raster_val: Bool. This is intended to offer methods that complement those used in Circuitscape.jl and Omniscape.jl. In this context, weights (the values in weight_raster) are in units of electrical resistance.  If false, the weight between two nodes with resistances R1 and R2 is  calculated by converting resistance to conductances, taking the average, then  taking the inverse of the result to convert back to resistance:  1 / ((1/R1 + 1/R2) / 2). connect_using_avg_weights = false correspondes to  the default settings in Circuitscape. Defaults to true', in which case the  simple average of the weights (adjusted for distance in the case of diagonal  neighbors) inweight_raster` are used.\n\ncombine: In the case that there are multiple edges defined for a single pair of vertices, how should the weight be chosen? Defaults to min. See the docs for SparseArrays.sparse() for more information.\n\n\n\n\n\n","category":"function"},{"location":"userguide/#Weighted-Graphs","page":"User Guide","title":"Weighted Graphs","text":"","category":"section"},{"location":"userguide/","page":"User Guide","title":"User Guide","text":"weightedrastergraph\nmake_weighted_raster_graph","category":"page"},{"location":"userguide/#SpatialGraphs.weightedrastergraph","page":"User Guide","title":"SpatialGraphs.weightedrastergraph","text":"weightedrastergraph(\n    weight_raster::Raster;\n    directed::Bool = false,\n    condition_raster::Raster = weight_raster,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false,\n    connect_using_avg_weights::Bool = true\n)\n\nConstruct a WeightedRasterGraph or WeightedRasterDiGraph (if  directed = true) from a raster dataset. The weight raster,  denotes the edge weight correponding to each vertex. Since edges are between rather than on vertices, edge weights are calculated as the average of the weights for each vertex.\n\nParameters\n\nweight_raster: A Rasters.Raster contained values that, where applicable  based on other arguments, determines which pixels to connect and the edge  weights between pixels. Any pixel in weight_raster with a value not equal to  weight_raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). \n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition_raster: A raster with values that can be used to determine whether two neighboring pixels should be connected. For example, an elevation raster  can be used to create a hydologic flow graph.\n\ncondition: A function that compares the values in condition_raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations with a lower value in condition_raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use < for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (condition_raster.missingval) in condition_raster. Note that if using an inequality function (or any function where the result depends on argument  position), the graph must be directed. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in condition_raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\nconnect_using_avg_weights: Bool. This is intended to offer methods that complement those used in Circuitscape.jl and Omniscape.jl. In this context, weights (the values in weight_raster) are in units of electrical resistance.  If false, the weight between two nodes with resistances R1 and R2 is  calculated by converting resistance to conductances, taking the average, then  taking the inverse of the result to convert back to resistance:  1 / ((1/R1 + 1/R2) / 2). connect_using_avg_weights = false corresponds to  the default settings in Circuitscape. Defaults to true, in which case the  simple average (adjusted for distance in the case of diagonal  neighbors) of the weights  in weight_raster is used.\n\n\n\n\n\n","category":"function"},{"location":"userguide/#SpatialGraphs.make_weighted_raster_graph","page":"User Guide","title":"SpatialGraphs.make_weighted_raster_graph","text":"make_weighted_raster_graph(\n    weight_raster::Raster,\n    vertex_raster::Raster;\n    directed::Bool = false,\n    condition_raster::Raster = weight_raster,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false,\n    connect_using_avg_weights::Bool = true,\n    combine::Function = min\n)\n\nConstruct a SimpleWeightedGraph or SimpleWeightedDiGraph (if  directed = true) based on vertex and weight rasters. This function is useful  if you already have a custom vertex raster and don't want SpatialGraphs.jl to make one for you. The vertex raster denotes the spatial locations of each vertex in the graph, and the weight raster denotes the edge weight correponding to each vertex. Since edges are between rather than on vertices, edge weights are calculated as the average of the weights for each vertex being connected.\n\nParameters\n\nweight_raster: A Rasters.Raster containing values that, where applicable  based on other arguments, determine which pixels to connect and the edge  weights between pixels. Any pixel in weight_raster with a value not equal to  weight_raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). \n\nvertex_raster: A Rasters.Raster with integer values ranging from 1:n,  where n is the number of unique vertices in the graph. \n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition: A function that compares the values in raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations that have a lower value in raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use > for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (raster.missingval). Note that if using an inequality function (or any function where the result depends on argument position), directed  should be set to true. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\nconnect_using_avg_weights: Bool. This is intended to offer methods that complement those used in Circuitscape.jl and Omniscape.jl. In this context, weights (the values in weight_raster) are in units of electrical resistance.  If false, the weight between two nodes with resistances R1 and R2 is  calculated by converting resistance to conductances, taking the average, then  taking the inverse of the result to convert back to resistance:  1 / ((1/R1 + 1/R2) / 2). connect_using_avg_weights = false correspondes to  the default settings in Circuitscape. Defaults to true, in which case the  simple average (adjusted for distance in the case of diagonal  neighbors) of the weights  in weight_raster is used.\n\ncombine: In the case that there are multiple edges defined for a single pair of vertices, how should the weight be chosen? Defaults to min. See the docs for SparseArrays.sparse() for more information.\n\n\n\n\n\n","category":"function"},{"location":"#SpatialGraphs.jl","page":"About","title":"SpatialGraphs.jl","text":"","category":"section"},{"location":"","page":"About","title":"About","text":"SpatialGraphs.jl introduces the AbstractSpatialGraph. AbstractSpatialGraphs  are a subtype of Graphs.AbstractGraph, and can be weighted or directed. SpatialGraphs.jl is useful for turning spatial data into graphs. This can be useful for landscape connectivity analysis, hydrology, and other spatial network processes. AbstractSpatialGraphs are AbstractGraphs, so methods from  Graphs.jl work right out of the box. Go to Graph Types  for more details on the graph types implemented in this package.","category":"page"}]
}
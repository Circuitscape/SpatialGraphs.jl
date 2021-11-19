var documenterSearchIndex = {"docs":
[{"location":"graphtypes/#graph_types","page":"Graph Types","title":"Graph Types in SpatialGraphs.jl","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"At this time, only raster-based graph types have been developed (and  vertex locations are stored in a Rasters.Raster), but there are plans to  eventually implement graph types for vector data as well.","category":"page"},{"location":"graphtypes/#Abstract-Types","page":"Graph Types","title":"Abstract Types","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"The overarching type in SpatialGraphs.jl is the AbstractSpatialGraph. All AbstractSpatialGraph subtypes contain a field called graph, which is itself an AbstractGraph.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"AbstractSpatialGraph","category":"page"},{"location":"graphtypes/#SpatialGraphs.AbstractSpatialGraph","page":"Graph Types","title":"SpatialGraphs.AbstractSpatialGraph","text":"AbstractSpatialGraph{T}\n\nAn abstract type representing a spatially referenced graph.\n\nAn AbstractSpatialGraph must contain the following field:\n\ngraph::AbstractGraph\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"The AbstractRasterGraph type is a subtype of AbstractSpatialGraph. All AbstractRasterGraph subtypes contain a field called vertex_raster, which is  a Rasters.Raster that describes the spatial locations for the vertices in  the graph. If your graph has 10 vertices, then the corresponding vertex_raster  must have 10 unique values, starting at 1 (which corresponds to the first vertex in the graph) and going up to 10 (which corresponds to the tenth vertex in the  graph). To find the location of a given graph vertex, n, you simply identify  the pixel(s) with a value of n. For pixels/elements in vertex_raster for which there shouldn't be a graph vertex, use a value of 0.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"AbstractRasterGraph","category":"page"},{"location":"graphtypes/#SpatialGraphs.AbstractRasterGraph","page":"Graph Types","title":"SpatialGraphs.AbstractRasterGraph","text":"AbstractRasterGraph{T}\n\nAn abstract type representing a spatially referenced graph, with graph vertices corresponding to pixels in a raster.\n\nAn AbstractRasterGraph must contain the following fields:\n\ngraph::AbstractGraph\nvertex_raster::Raster\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#AbstractRasterGraph-Subtypes","page":"Graph Types","title":"AbstractRasterGraph Subtypes","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"SpatialGraphs.jl has raster graph types for undirected, directed, unweighted, and weighted graphs, each detailed below.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"RasterGraph\nRasterDiGraph\nWeightedRasterGraph\nWeightedRasterDiGraph","category":"page"},{"location":"graphtypes/#SpatialGraphs.RasterGraph","page":"Graph Types","title":"SpatialGraphs.RasterGraph","text":"RasterGraph{T}\n\nA composite type for a spatially referenced graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#SpatialGraphs.RasterDiGraph","page":"Graph Types","title":"SpatialGraphs.RasterDiGraph","text":"RasterDiGraph{T}\n\nA composite type for a spatially referenced directed graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#SpatialGraphs.WeightedRasterGraph","page":"Graph Types","title":"SpatialGraphs.WeightedRasterGraph","text":"WeightedRasterGraph{T}\n\nA composite type for a spatially referenced weighted graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"graphtypes/#SpatialGraphs.WeightedRasterDiGraph","page":"Graph Types","title":"SpatialGraphs.WeightedRasterDiGraph","text":"WeightedRasterDiGraph{T}\n\nA composite type for a spatially referenced, weighted, directed graph. Vertices are spatially referenced based on a raster.\n\n\n\n\n\n","category":"type"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Habitat-Corridor-Mapping","page":"Examples","title":"Habitat Corridor Mapping","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this example, we'll map a habitat corridor in central Maryland. We'll use SpatialGraphs.jl to convert the landscape into a graph, then calculate cost  distances, and finally convert the result back into a raster to display a least cost corridor in space.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"First, install the necessary packages and import them:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Pkg\nPkg.add([\"Rasters\", \"Plots\", \"Graphs\", \"SimpleWeightedGraphs\", \"SpatialGraphs\"])\nusing SpatialGraphs, Rasters, Plots, Graphs, SimpleWeightedGraphs","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using Pkg\nPkg.add([\"Rasters\", \"Graphs\", \"SimpleWeightedGraphs\"])\nusing SpatialGraphs, Rasters, Graphs, SimpleWeightedGraphs\n\nurl_base = \"https://raw.githubusercontent.com/Circuitscape/datasets/main/\"\ndownload(string(url_base, \"data/nlcd_2016_frederick_md.tif\"),\n         \"nlcd_2016_frederick_md.tif\")\n\nland_cover = Raster(\"nlcd_2016_frederick_md.tif\")\nnothing","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now, let's download the land cover map that we'll use as the basis for defining the weights of our graph, and plot it. Graph weights correspond to the cost of  an animal moving from one graph vertex to another. This is also commonly  referred to as \"resistance\" in the context of connectivity modeling.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"url_base = \"https://raw.githubusercontent.com/Circuitscape/datasets/main/\"\ndownload(string(url_base, \"nlcd_2016_frederick_md.tif\"),\n         \"nlcd_2016_frederick_md.tif\")\nvalues = [11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95]\npalette = [\"#476BA0\", \"#DDC9C9\", \"#D89382\", \"#ED0000\", \"#AA0000\",\n           \"#b2b2b2\", \"#68AA63\", \"#1C6330\", \"#B5C98E\", \"#CCBA7C\",\n           \"#E2E2C1\", \"#DBD83D\", \"#AA7028\", \"#BAD8EA\", \"#70A3BA\"]\n\nplot(Raster(\"nlcd_2016_frederick_md.tif\"), title = \"Land Cover Type\", \n     xlabel = \"Easting\", ylabel = \"Northing\", \n     seriescolor = cgrad(palette, (values .- 12) ./ 84,\n     categorical = true), size = (700, 640))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<img src='../figs/mdlc.png' width=500><br>","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Load the land cover data as a Raster, read it into memory, then convert it to  Float64.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"landcover = convert.(\n    Float64,\n    read(\n        Raster(\n            \"nlcd_2016_frederick_md.tif\",\n            missingval = -9999.0\n        )\n    )\n)\n","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now, we need to convert land cover into resistance by reclassifying the values in the land cover raster to their corresponding costs/weights.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"# Copy the landcover object to use as a place holder for resistance\nresistance = deepcopy(landcover)\n\n# Define how values should be reclassified. Original values are on the left,\n# and new values are on the right.\nreclass_table = [\n    11.\t100; # Water\n    21\t500; # Developed, open space\n    22\t1000; # Developed, low intensity\n    23\t-9999; # Developed, medium intensity\n    24\t-9999; # Developed, high intensity\n    31\t100; # Barren land\n    41\t1; # Deciduous forest\n    42\t1; # Evergreen forest\n    43\t1; # Mixed forest\n    52\t20; # Shrub/scrub\n    71\t30; # Grassland/herbaceous\n    81\t200; # Pasture/hay\n    82\t300; # Cultivated crops\n    90\t20; # Woody wetlands\n    95\t30; # Emergent herbaceous wetlands\n]\n\n# Redefine the values in resistance based on the corresponding land cover\n# values, using reclass_table \nfor i in 1:(size(reclass_table)[1])\n    resistance[landcover .== reclass_table[i, 1]] .= reclass_table[i, 2]\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now that we have a resistance raster, we can start using SpatialGraphs.jl!","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Convert the resistance raster into a WeightedRasterGraph:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"wrg = weightedrastergraph(resistance)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now that we have a WeightedRasterGraph, we can make use of the functions in Graphs.jl.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Let's calculate cost distances from two different vertices. Vertex IDs,  the second argument to dijkstra_shortest_paths, correspond to the values in  wrg.vertex_raster.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"First, we'll use dijkstra_shortest_paths to calculate to total cost distance from each vertex of interest to all other vertices in the graph.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"cd_1 = dijkstra_shortest_paths(wrg, 1)\ncd_2 = dijkstra_shortest_paths(wrg, 348021)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now, create placeholder rasters to store the cost distance values for each  vertex, and populate them using the cost distances calculated above. Finally, we'll sum the rasters to enable us to delineate a corridor.  The value at each pixel corresponds to the total cost distance of the least cost path between vertices 1 and 348021 that intersects that pixel.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"# Create rasters\ncd_1_raster = Raster( # Start with an empty raster that will be populated\n    fill(float(resistance.missingval), size(wrg.vertex_raster)),\n    dims(resistance),\n    missingval=resistance.missingval\n)\ncd_2_raster = Raster(\n    fill(float(wrg.vertex_raster.missingval), size(wrg.vertex_raster)),\n    dims(resistance),\n    missingval=resistance.missingval\n)\n\n# Populate with cost distances\ncd_1_raster[wrg.vertex_raster .!= wrg.vertex_raster.missingval] .= cd_1.dists[\n    wrg.vertex_raster[wrg.vertex_raster.!= wrg.vertex_raster.missingval]\n]\ncd_2_raster[wrg.vertex_raster .!= wrg.vertex_raster.missingval] .= cd_2.dists[\n    wrg.vertex_raster[wrg.vertex_raster.!= wrg.vertex_raster.missingval]\n]\n\n# Sum the rasters to enalbe corridor mapping\ncwd_sum = cd_1_raster + cd_2_raster","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Setting some land cover classes to NoData resulted in some vertices being completely isolated from the vertices from which we calculated cost distances above. In these cases, dijkstra_shortest_paths returns a value of Inf for the cost distance. Let's get rid of those by setting them to the NoData value of our raster.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"cwd_sum[cwd_sum .== Inf] .= cwd_sum.missingval","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Finally, let's plot the result! We'll use a threshold to set pixels above a certain value to NoData, so we will only retain the pixels within our \"corridor\".","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"corridor = copy(cwd_sum)\ncorridor[cwd_sum .>= (minimum(cwd_sum[cwd_sum .>= 0]) + 250)] .= -9999\n\nplot(\n    cwd_sum,\n    title = \"Cost Distance Corridor\",\n    xlabel = \"Easting\", ylabel = \"Northing\",\n    colorbar_title = \" \\nCost Distance Sum\",\n    size = (650, 475)\n)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<img src='../figs/corridor_plot.png' width=500><br>","category":"page"},{"location":"userguide/#User-Guide","page":"User Guide","title":"User Guide","text":"","category":"section"},{"location":"userguide/#Building-Graphs-from-Rasters","page":"User Guide","title":"Building Graphs from Rasters","text":"","category":"section"},{"location":"userguide/","page":"User Guide","title":"User Guide","text":"SpatialGraphs.jl offers several functions for constructing graphs from raster data.","category":"page"},{"location":"userguide/#Simple-Graphs","page":"User Guide","title":"Simple Graphs","text":"","category":"section"},{"location":"userguide/","page":"User Guide","title":"User Guide","text":"rastergraph\nmake_simple_raster_graph","category":"page"},{"location":"userguide/#SpatialGraphs.rastergraph","page":"User Guide","title":"SpatialGraphs.rastergraph","text":"RasterGraph(\n    raster::Raster;\n    directed::Bool = true,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false\n)\n\nConstruct a RasterGraph or RasterDiGraph (if  directed = true) from a raster dataset.\n\nParameters\n\nraster: A Rasters.Raster on which to base the graph. Any pixel in raster  with a value not equal to raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). The values in the raster can also be used to determine which vertices to connect. See condition below for more information.\n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition: A function that compares the values in raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations that have a lower value in raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use > for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (raster.missingval). Note that if using an inequality function (or any function where the result depends on argument position), directed  should be set to true. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\n\n\n\n\n","category":"function"},{"location":"userguide/#SpatialGraphs.make_simple_raster_graph","page":"User Guide","title":"SpatialGraphs.make_simple_raster_graph","text":"make_simple_raster_graph(\n    raster::Raster,\n    vertex_raster::Raster;\n    directed::Bool = false,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false,\n    combine::Function = min\n)\n\nConstruct a SimpleGraph or SimpleDiGraph (if directed = true) based on  two rasters (raster and vertex_raster). This function is useful  if you already have a custom vertex raster and don't want SpatialGraphs.jl to make one for you. The vertex raster denotes the spatial locations of each vertex in the graph, and raster is used to construct the graph and determine which vertices to connect.\n\nParameters\n\nraster: A Rasters.Raster on which to base the graph. Any pixel in raster  with a value not equal to raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). The values in the raster can also be used to determine which vertices to connect. See condition below for more information.\n\nvertex_raster: A Rasters.Raster with integer values ranging from 1:n,  where n is the number of unique vertices in the graph. \n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition: A function that compares the values in condition_raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations with a lower value in condition_raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use < for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (condition_raster.missingval) in condition_raster. Note that if using an inequality function (or any function where the result depends on argument  position), the graph must be directed. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in condition_raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\nconnect_using_avg_raster_val: Bool. This is intended to offer methods that complement those used in Circuitscape.jl and Omniscape.jl. In this context, weights (the values in weight_raster) are in units of electrical resistance.  If false, the weight between two nodes with resistances R1 and R2 is  calculated by converting resistance to conductances, taking the average, then  taking the inverse of the result to convert back to resistance:  1 / ((1/R1 + 1/R2) / 2). connect_using_avg_weights = false correspondes to  the default settings in Circuitscape. Defaults to true', in which case the  simple average of the weights (adjusted for distance in the case of diagonal  neighbors) inweight_raster` are used.\n\ncombine: In the case that there are multiple edges defined for a single pair of vertices, how should the weight be chosen? Defaults to min. See the docs for SparseArrays.sparse() for more information.\n\n\n\n\n\n","category":"function"},{"location":"userguide/#Weighted-Graphs","page":"User Guide","title":"Weighted Graphs","text":"","category":"section"},{"location":"userguide/","page":"User Guide","title":"User Guide","text":"weightedrastergraph\nmake_weighted_raster_graph","category":"page"},{"location":"userguide/#SpatialGraphs.weightedrastergraph","page":"User Guide","title":"SpatialGraphs.weightedrastergraph","text":"weightedrastergraph(\n    weight_raster::Raster;\n    directed::Bool = false,\n    condition_raster::Raster = weight_raster,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false,\n    connect_using_avg_weights::Bool = true\n)\n\nConstruct a WeightedRasterGraph or WeightedRasterDiGraph (if  directed = true) from a raster dataset. The weight raster,  denotes the edge weight correponding to each vertex. Since edges are between rather than on vertices, edge weights are calculated as the average of the weights for each vertex.\n\nParameters\n\nweight_raster: A Rasters.Raster contained values that, where applicable  based on other arguments, determines which pixels to connect and the edge  weights between pixels. Any pixel in weight_raster with a value not equal to  weight_raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). \n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition_raster: A raster with values that can be used to determine whether two neighboring pixels should be connected. For example, an elevation raster  can be used to create a hydologic flow graph.\n\ncondition: A function that compares the values in condition_raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations with a lower value in condition_raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use < for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (condition_raster.missingval) in condition_raster. Note that if using an inequality function (or any function where the result depends on argument  position), the graph must be directed. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in condition_raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\nconnect_using_avg_weights: Bool. This is intended to offer methods that complement those used in Circuitscape.jl and Omniscape.jl. In this context, weights (the values in weight_raster) are in units of electrical resistance.  If false, the weight between two nodes with resistances R1 and R2 is  calculated by converting resistance to conductances, taking the average, then  taking the inverse of the result to convert back to resistance:  1 / ((1/R1 + 1/R2) / 2). connect_using_avg_weights = false corresponds to  the default settings in Circuitscape. Defaults to true, in which case the  simple average (adjusted for distance in the case of diagonal  neighbors) of the weights  in weight_raster is used.\n\n\n\n\n\n","category":"function"},{"location":"userguide/#SpatialGraphs.make_weighted_raster_graph","page":"User Guide","title":"SpatialGraphs.make_weighted_raster_graph","text":"make_weighted_raster_graph(\n    weight_raster::Raster,\n    vertex_raster::Raster;\n    directed::Bool = false,\n    condition_raster::Raster = weight_raster,\n    condition::Function = is_data,\n    cardinal_neighbors_only::Bool = false,\n    connect_using_avg_weights::Bool = true,\n    combine::Function = min\n)\n\nConstruct a SimpleWeightedGraph or SimpleWeightedDiGraph (if  directed = true) based on vertex and weight rasters. This function is useful  if you already have a custom vertex raster and don't want SpatialGraphs.jl to make one for you. The vertex raster denotes the spatial locations of each vertex in the graph, and the weight raster denotes the edge weight correponding to each vertex. Since edges are between rather than on vertices, edge weights are calculated as the average of the weights for each vertex being connected.\n\nParameters\n\nweight_raster: A Rasters.Raster containing values that, where applicable  based on other arguments, determine which pixels to connect and the edge  weights between pixels. Any pixel in weight_raster with a value not equal to  weight_raster.missingval will be assigned a vertex in the graph (corresponding to its centroid). \n\nvertex_raster: A Rasters.Raster with integer values ranging from 1:n,  where n is the number of unique vertices in the graph. \n\nArguments\n\ndirected: A Bool determining whether the graph should be directed.\n\ncondition: A function that compares the values in raster for two neighbors to determine if those neighbors should be connected. The function must compare two values and return either true or false. Useful functions to use  here include <, <=, ==, etc. The first argument to condition corresponds to the source vertex, and the second argument corresponds to the destination  vertex. So, if you only want to connect sources to destinations that have a lower value in raster (e.g. in the case of developing a hydrologic flow  graph based on elevation), then you would use > for condition. Defaults to  is_data, which results in neighbors being connected as long as they are not  NoData (raster.missingval). Note that if using an inequality function (or any function where the result depends on argument position), directed  should be set to true. For undirected graphs, you can use either is_data or ==, or any other custom function where argument position doesn't matter, e.g. a function that determines whether the values in raster are within a certain distance of each other.\n\ncardinal_neighbors_only: A Bool stating whether only cardinal neighbors should be connected. By default, both cardinal and diagonal neighbors are connected. Note that when determining weights between diagonal neighbors, the increased distance between them (as compared to the distance between cardinal neighbors) is accounted for.\n\nconnect_using_avg_weights: Bool. This is intended to offer methods that complement those used in Circuitscape.jl and Omniscape.jl. In this context, weights (the values in weight_raster) are in units of electrical resistance.  If false, the weight between two nodes with resistances R1 and R2 is  calculated by converting resistance to conductances, taking the average, then  taking the inverse of the result to convert back to resistance:  1 / ((1/R1 + 1/R2) / 2). connect_using_avg_weights = false correspondes to  the default settings in Circuitscape. Defaults to true, in which case the  simple average (adjusted for distance in the case of diagonal  neighbors) of the weights  in weight_raster is used.\n\ncombine: In the case that there are multiple edges defined for a single pair of vertices, how should the weight be chosen? Defaults to min. See the docs for SparseArrays.sparse() for more information.\n\n\n\n\n\n","category":"function"},{"location":"#SpatialGraphs.jl","page":"About","title":"SpatialGraphs.jl","text":"","category":"section"},{"location":"","page":"About","title":"About","text":"SpatialGraphs.jl introduces the AbstractSpatialGraph. AbstractSpatialGraphs  are a subtype of Graphs.AbstractGraph, and can be weighted or directed. SpatialGraphs.jl is useful for turning spatial data into graphs. This can be useful for landscape connectivity analysis, hydrology, and other spatial network processes. AbstractSpatialGraphs are AbstractGraphs, so methods from  Graphs.jl work right out of the box. Go to Graph Types  for more details on the graph types implemented in this package.","category":"page"},{"location":"#Table-of-Contents","page":"About","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"About","title":"About","text":"Pages = [\"graphtypes.md\",\"userguide.md\", \"examples.md\"]\nDepth = 2","category":"page"}]
}

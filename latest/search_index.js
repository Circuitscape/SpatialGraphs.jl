var documenterSearchIndex = {"docs":
[{"location":"graphtypes/#Graph-Types-in-SpatialGraphs.jl","page":"Graph Types","title":"Graph Types in SpatialGraphs.jl","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"At this time, only raster-based graph types have been developed (and  vertex locations are stored in a GeoData.GeoArray), but there are plans to  eventually implement graph types for vector data as well.","category":"page"},{"location":"graphtypes/#Abstract-Types","page":"Graph Types","title":"Abstract Types","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"The overarching type in SpatialGraphs.jl is the AbstractSpatialGraph. All AbstractSpatialGraph subtypes contain a field called graph, which is itself an AbstractGraph.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"AbstractSpatialGraph","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"The AbstractRasterGraph type is as subtype of AbstractSpatialGraph. All AbstractRasterGraph subtypes contain a field called vertex_raster, which is  a GeoData.GeoArray that describes the spatial locations for the vertices in  the graph. If your graph has 10 vertices, then the corresponding GeoArray  must have 10 unique values, starting at 1 (which corresponds to the first vertex in the graph) and going up to 10 (which corresponds to the tenth vertex in the  graph). To find the location of a given graph vertex, n, you simply identify  the pixel(s) with a value of n.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"AbstractRasterGraph","category":"page"},{"location":"graphtypes/#AbstractRasterGraph-Subtypes","page":"Graph Types","title":"AbstractRasterGraph Subtypes","text":"","category":"section"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"SpatialGraphs.jl has raster graph types for undirected, directed, unweighted, and weighted graphsm each detailed below.","category":"page"},{"location":"graphtypes/","page":"Graph Types","title":"Graph Types","text":"SimpleRasterGraph\nSimpleRasterDiGraph\nWeightedRasterGraph\nWeightedRasterDiGraph","category":"page"},{"location":"#SpatialGraphs.jl","page":"About","title":"SpatialGraphs.jl","text":"","category":"section"},{"location":"","page":"About","title":"About","text":"SpatialGraphs.jl introduces the AbstractSpatialGraph. AbstractSpatialGraphs  are a subtype of LightGraphs.AbstractGraph, and can be weighted or directed. SpatialGraphs.jl is useful for turning spatial data into graphs. This can be useful for landscape connectivity analysis, hydrology, and other spatial network processes. AbstractSpatialGraphs are AbstractGraphs, methods from  LightGraphs.jl and SimpleWeightedGraphs.jl work right out of the box.  Go to Graph Types for more details on the graph types implemented in  this package.","category":"page"}]
}

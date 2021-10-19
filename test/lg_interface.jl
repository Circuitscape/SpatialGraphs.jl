using SpatialGraphs, Test, GeoData, Graphs, SimpleWeightedGraphs
## This script tests methods for the LightGraph interface that aren't already 
## used in other tests

A_array = Array{Float32}(undef, (3, 4, 1))
A_array[:,:,:] = [1, 3, 2, 0.5, 10, 8, 5, -9999, 3, 1, 2, 6]
x = X(1:4)
y = Y(1:3)
band = Band(1:1)

weight_raster = GeoArray(A_array, (y, x, band), missingval = -9999)
rasgraph = weightedrastergraph(weight_raster)

@test nv(rasgraph) == maximum(rasgraph.vertex_raster)
num_edge = ne(rasgraph)
@test vertices(rasgraph) == 1:maximum(rasgraph.vertex_raster)
@test eltype(rasgraph) == Int # Should always be Int64 (or int32 on 32-bit)
@test edgetype(rasgraph) == SimpleWeightedEdge{eltype(rasgraph), eltype(A_array)}

for i in 1:maximum(rasgraph.vertex_raster)
    @test has_vertex(rasgraph, i)
end

# In undirected graph incoming neighbors should equal outgoing
@test inneighbors(rasgraph, 1) == outneighbors(rasgraph, 1)

empty = zero(rasgraph)

add_vertices!(empty, maximum(rasgraph.vertex_raster))
add_edge!(empty, 1, 4, Float32(0.5))

@test nv(empty) == maximum(rasgraph.vertex_raster)
@test has_edge(empty, 1, 4)
@test get_weight(empty, 1, 4) == Float32(0.5)
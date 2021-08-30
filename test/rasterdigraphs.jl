using ArchGDAL, GeoData, LightGraphs, SimpleWeightedGraphs, SpatialGraphs
A_array = Array{Float64}(undef, (3, 4, 1))
A_array[:,:,:] = [1, 3, 2, 0.5, 10, 8, 5, -9999, 3, 1, 2, 6]

condition_array = Array{Float64}(undef, (3, 4, 1))
condition_array[:,:,:] = [1, 3, 5, 2, 4, 8, 5, -9999, 2, 3, 6, 7]

x = X(1:4)
y = Y(1:3)
band = Band(1:1)

weight_raster = GeoArray(A_array, (y, x, band), missingval = -9999)
condition_raster = GeoArray(condition_array, (y, x, band), missingval = -9999)

compare = <
rasgraph = weightedrastergraph(
    weight_raster,
    directed = true,
    condition_raster = condition_raster,
    condition = compare
)

# no vertices in NoData pixels?
@test (rasgraph.vertex_raster .== 0) == 
        ((weight_raster .== weight_raster.missingval) .| 
         isnan.(weight_raster))

# Is the number of of vertices correct, and is the range of values correct?
@test sort(collect(rasgraph.vertex_raster[rasgraph.vertex_raster .!= 0])) == 
        collect(1:sum(
            (weight_raster .!= weight_raster.missingval) .& 
                (!).(isnan.(weight_raster))
        ))

graph_edges = collect(edges(rasgraph))

# Test that the edges are correct and have proper weights
for i in 1:length(graph_edges)
    source_i = src(graph_edges[i])
    dest_i = dst(graph_edges[i])
    weight_i = graph_edges[i].weight

    source_coords = findall(rasgraph.vertex_raster .== source_i)[1]
    dest_coords = findall(rasgraph.vertex_raster .== dest_i)[1]

    # Check that condition is met
    @test compare(
        condition_raster[source_coords],
        condition_raster[dest_coords]
    )

    # Test that source row is within 1 step of dest row
    row_diff = abs(source_coords[1] - dest_coords[1])
    @test row_diff <= 1

    # Test that source column is within 1 step of dest row
    col_diff = abs(source_coords[2] - dest_coords[2])
    @test col_diff <= 1

    # Test that the weight is what it should be (assumes connect_using_avg_weights = true in graph construction)
    if (row_diff == 1 && col_diff == 1) # get diagonal average
        @test weight_i == 
            SpatialGraphs.res_diagonal_avg(weight_raster[source_coords], 
                                           weight_raster[dest_coords])
    else
        @test weight_i == 
            SpatialGraphs.res_cardinal_avg(weight_raster[source_coords],
                                           weight_raster[dest_coords])
    end
end

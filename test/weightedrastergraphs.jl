using GeoData, LightGraphs, SimpleWeightedGraphs, SpatialGraphs, Test
A_array = Array{Float64}(undef, (3, 4, 1))
A_array[:,:,:] = [1, 3, 2, 0.5, 10, 8, 5, -9999, 3, 1, 2, 6]

x = X(1:4)
y = Y(1:3)
band = Band(1:1)

weight_raster = GeoArray(A_array, (y, x, band), missingval = -9999)
weight_raster2d = GeoArray(A_array[:, :, 1], (y, x), missingval = -9999)
rasgraph2d = weightedrastergraph(weight_raster2d)
rasgraph = weightedrastergraph(weight_raster)

# Test that graphs are the same regardless of whether weight_raster has Band dim
@test rasgraph2d.graph == rasgraph.graph

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

# check that no edges exist that shouldn't
vertex_raster = rasgraph.vertex_raster
for src_row in 1:size(vertex_raster)[1]
    for src_col in 1:size(vertex_raster)[2]
        for dst_row in 1:size(vertex_raster)[1]
            for dst_col in 1:size(vertex_raster)[2]
                if ((vertex_raster[src_row, src_col] == 0) | 
                        (vertex_raster[dst_row, dst_col] == 0))
                    continue
                elseif ((abs(src_row - dst_row) <= 1) &
                        (abs(src_col - dst_col) <= 1)) & (
                            (abs(src_row - dst_row) > 0) |
                            (abs(src_col - dst_col) > 0))
                    @test has_edge(
                        rasgraph, 
                        vertex_raster[src_row, src_col],
                        vertex_raster[dst_row, dst_col]
                    ) == true
                else
                    @test (has_edge(
                        rasgraph, 
                        vertex_raster[src_row, src_col],
                        vertex_raster[dst_row, dst_col]
                    ) == false)
                end
            end
        end
    end
end

@test is_directed(rasgraph) == false
@test zero(rasgraph).vertex_raster == rasgraph.vertex_raster
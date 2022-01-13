using Rasters, SpatialGraphs, Test, Graphs, SimpleWeightedGraphs

resistance = [
    1. 2.5 6 0.2;
    5 6 12 1;
    0.5 3 7 4;
    11 0.5 9 -9999
]

res = Raster(
    reshape(resistance, (size(resistance)..., 1)),
    missingval=-9999,
    dims=(X(1:4), Y(1:4), Band(1:1))
)

patches = [
    1 1 0 3;
    0 0 0 3;
    0 0 0 0;
    0 2 2 0
]

pat = Raster(
    reshape(patches, (size(patches)..., 1)),
    missingval=-9999,
    dims=(X(1:4), Y(1:4), Band(1:1))
)

vertex_raster = make_vertex_raster(res, pat)

# Check that patches have proper vertex ID
for i in 1:maximum(pat)
   @test all(vertex_raster[pat .== i] .== i)
end

@test sort(unique(vertex_raster[vertex_raster .> maximum(pat)])) ==
    (maximum(pat) + 1):(maximum(pat) + length(res[(pat .== 0) .& (res .!= res.missingval)]))

## Make sure that a graph can be successfully created from a the vertex raster
g = WeightedRasterGraph(make_weighted_raster_graph(res, vertex_raster), vertex_raster)

# Test that a few random weights are correct, assumes that defaults were used
# in make_weighted_raster_graph above, particularly that combine = min and
# connect_using_avg_weights = true. THESE ARE HARD CODED. Changes to res or pat
# above will result in failed tests.
@test get_weight(g.graph, 1, 4) == 3 # minimum of the mean resistances between both vertices labeled 1 and vertex labeled 4
@test get_weight(g.graph, 7, 11) == SpatialGraphs.res_diagonal_avg(res[vertex_raster .== 7], res[vertex_raster .== 11])[1]
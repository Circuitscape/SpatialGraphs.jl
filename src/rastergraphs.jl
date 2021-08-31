"""
    make_vertex_raster(A::GeoArray)

Constuct a vertex raster (a raster where the value of each pixel corresponds
to its ID in a graph, and 0s correspond to NoData). Returns a GeoArray. This 
function is recommended for internal use only.

## Parameters
`A`: The GeoArray from which a graph will be built, which is used as the
reference for building the vertex raster. Pixels with NoData (`A.missingval`)
are skipped (no vertex is assigned). Pixels with NoData will get a value of 0 in
the vertex raster.
"""
function make_vertex_raster(A::GeoData.GeoArray)
    # Make an array of unique node identifiers
    nodemap = zeros(Int64, size(A.data))
    is_node = (A.data .!= A.missingval) .&
        ((!).(isnan.(A.data)))

    nodemap[is_node] = 1:sum(is_node)

    nodemap = GeoData.GeoArray(nodemap, dims(A))

    nodemap
end


"""
    weightedrastergraph(
        weight_raster::GeoArray;
        directed::Bool = false,
        condition_raster::GeoArray = weight_raster,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        connect_using_avg_weights::Bool = true
    )

Construct a `WeightedRasterGraph` or `WeightedRasterDiGraph` (if 
`directed = true`) from a raster dataset. The weight raster,  denotes
the edge weight correponding to each vertex. Since edges are between rather than
on vertices, edge weights are calculated as the average of the weights for each
vertex.

## Parameters
`weight_raster`: A GeoData.GeoArray contained values that, where applicable 
based on other arguments, determines which pixels to connect and the edge 
weights between pixels. Any pixel in `weight_raster` with a value not equal to 
`weight_raster.missingval` will be assigned a vertex in the graph (corresponding
to its centroid). 

## Arguments
`directed`: A `Bool` determining whether the graph should be directed.

`condition_raster`: A raster with values that can be used to determine whether
two neighboring pixels should be connected. For example, an elevation raster 
can be used to create a hydologic flow graph.

`condition`: A function that compares the values in `condition_raster` for two
neighbors to determine if those neighbors should be connected. The function must
compare two values and return either `true` or `false`. Useful functions to use 
here include `<`, `<=`, `==`, etc. The first argument to `condition` corresponds
to the source vertex, and the second argument corresponds to the destination 
vertex. So, if you only want to connect sources to destinations with a lower
value in `condition_raster` (e.g. in the case of developing a hydrologic flow 
graph based on elevation), then you would use `<` for `condition`. Defaults to 
`is_data`, which results in neighbors being connected as long as they are not 
NoData (`condition_raster.missingval`) in `condition_raster`. Note that if using
an inequality function (or any function where the result depends on argument 
position), the graph must be directed. For undirected graphs, you can use either
`is_data` or `==`, or any other custom function where argument position doesn't
matter, e.g. a function that determines whether the values in `condition_raster`
are within a certain distance of each other.

`cardinal_neighbors_only`: A `Bool` stating whether only cardinal neighbors
should be connected. By default, both cardinal _and_ diagonal neighbors are
connected. Note that when determining weights between diagonal neighbors, the
increased distance between them (as compared to the distance between cardinal
neighbors) is accounted for.

`connect_using_avg_weights`: `Bool`. This is intended to offer methods that
complement those used in Circuitscape.jl and Omniscape.jl. In this context,
weights (the values in `weight_raster`) are in units of electrical resistance. 
If `false`, the weight between two nodes with resistances R1 and R2 is 
calculated by converting resistance to conductances, taking the average, then 
taking the inverse of the result to convert back to resistance: 
`1 / ((1/R1 + 1/R2) / 2)`. `connect_using_avg_weights = false` correspondes to 
the default settings in Circuitscape. Defaults to `true`, in which case the 
simple average (adjusted for distance in the case of diagonal 
neighbors) of the weights  in `weight_raster` is used.

"""
function weightedrastergraph(
        weight_raster::GeoArray;
        directed::Bool = false,
        condition_raster::GeoArray = weight_raster,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        connect_using_avg_weights::Bool = true
    )
    v = make_vertex_raster(weight_raster)
    g = make_weighted_raster_graph(
        weight_raster,
        v,
        directed = directed,
        condition_raster = condition_raster,
        condition = condition,
        cardinal_neighbors_only = cardinal_neighbors_only,
        connect_using_avg_weights = connect_using_avg_weights
    )

    if directed
        sg = WeightedRasterDiGraph(g, v)
    else
        sg = WeightedRasterGraph(g, v)
    end

    return sg

end

"""
    simplerastergraph(
        raster::GeoArray;
        directed::Bool = true,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false
    )

Construct a `SimpleRasterGraph` or `SimpleRasterDiGraph` (if 
`directed = true`) from a raster dataset.

## Parameters
`raster`: A GeoData.GeoArray on which to base the graph. Any pixel in `raster` 
with a value not equal to `raster.missingval` will be assigned a vertex
in the graph (corresponding to its centroid). The values in the raster can also
be used to determine which vertices to connect. See `condition` below for more
information.

## Arguments
`directed`: A `Bool` determining whether the graph should be directed.

`condition`: A function that compares the values in `raster` for two
neighbors to determine if those neighbors should be connected. The function must
compare two values and return either `true` or `false`. Useful functions to use 
here include `<`, `<=`, `==`, etc. The first argument to `condition` corresponds
to the source vertex, and the second argument corresponds to the destination 
vertex. So, if you only want to connect sources to destinations that have a
lower value in `raster` (e.g. in the case of developing a hydrologic flow 
graph based on elevation), then you would use `>` for `condition`. Defaults to 
`is_data`, which results in neighbors being connected as long as they are not 
NoData (`raster.missingval`). Note that if using an inequality function (or
any function where the result depends on argument position), `directed` 
should be set to `true`. For undirected graphs, you can use either
`is_data` or `==`, or any other custom function where argument position doesn't
matter, e.g. a function that determines whether the values in `raster`
are within a certain distance of each other.

`cardinal_neighbors_only`: A `Bool` stating whether only cardinal neighbors
should be connected. By default, both cardinal _and_ diagonal neighbors are
connected. Note that when determining weights between diagonal neighbors, the
increased distance between them (as compared to the distance between cardinal
neighbors) is accounted for.
"""
function simplerastergraph(
    raster::GeoArray;
    condition::Function = is_data,
    directed::Bool = true,
    cardinal_neighbors_only::Bool = false,
)
    v = make_vertex_raster(weight_raster)
    g = make_simple_raster_graph(
        raster,
        vertex_raster,
        directed = directed,
        condition = condition,
        cardinal_neighbors_only = cardinal_neighbors_only
    )

    if directed
        sg = SimpleRasterDiGraph(g, v)
    else
        sg = SimpleRasterGraph(g, v)
    end

    return sg
end

"""
    make_weighted_raster_graph(
        weight_raster::GeoArray,
        vertex_raster::GeoArray;
        directed::Bool = false,
        condition_raster::GeoArray = weight_raster,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        connect_using_avg_weights::Bool = true,
        combine::Function = min
    )

Construct a `SimpleWeightedGraph` or `SimpleWeightedDiGraph` (if 
`directed = true`) based on vertex and weight rasters. This function is useful 
if you already have a custom vertex raster and don't want SpatialGraphs.jl to
make one for you. The vertex raster denotes the spatial locations of each vertex
in the graph, and the weight raster denotes the edge weight correponding to each
vertex. Since edges are between rather than on vertices, edge weights are
calculated as the average of the weights for each vertex being connected.

## Parameters
`weight_raster`: A `GeoData.GeoArray` containing values that, where applicable 
based on other arguments, determine which pixels to connect and the edge 
weights between pixels. Any pixel in `weight_raster` with a value not equal to 
`weight_raster.missingval` will be assigned a vertex in the graph (corresponding
to its centroid). 

`vertex_raster`: A `GeoData.GeoArray` with integer values ranging from 1:n, 
where n is the number of unique vertices in the graph. 

## Arguments
`directed`: A `Bool` determining whether the graph should be directed.

`condition`: A function that compares the values in `raster` for two
neighbors to determine if those neighbors should be connected. The function must
compare two values and return either `true` or `false`. Useful functions to use 
here include `<`, `<=`, `==`, etc. The first argument to `condition` corresponds
to the source vertex, and the second argument corresponds to the destination 
vertex. So, if you only want to connect sources to destinations that have a
lower value in `raster` (e.g. in the case of developing a hydrologic flow 
graph based on elevation), then you would use `>` for `condition`. Defaults to 
`is_data`, which results in neighbors being connected as long as they are not 
NoData (`raster.missingval`). Note that if using an inequality function (or
any function where the result depends on argument position), `directed` 
should be set to `true`. For undirected graphs, you can use either
`is_data` or `==`, or any other custom function where argument position doesn't
matter, e.g. a function that determines whether the values in `raster`
are within a certain distance of each other.

`cardinal_neighbors_only`: A `Bool` stating whether only cardinal neighbors
should be connected. By default, both cardinal _and_ diagonal neighbors are
connected. Note that when determining weights between diagonal neighbors, the
increased distance between them (as compared to the distance between cardinal
neighbors) is accounted for.

`connect_using_avg_weights`: `Bool`. This is intended to offer methods that
complement those used in Circuitscape.jl and Omniscape.jl. In this context,
weights (the values in `weight_raster`) are in units of electrical resistance. 
If `false`, the weight between two nodes with resistances R1 and R2 is 
calculated by converting resistance to conductances, taking the average, then 
taking the inverse of the result to convert back to resistance: 
`1 / ((1/R1 + 1/R2) / 2)`. `connect_using_avg_weights = false` correspondes to 
the default settings in Circuitscape. Defaults to `true`, in which case the 
simple average (adjusted for distance in the case of diagonal 
neighbors) of the weights  in `weight_raster` is used.

`combine`: In the case that there are multiple edges defined for a single pair
of vertices, how should the weight be chosen? Defaults to `min`. See the docs
for `SparseArrays.sparse()` for more information.
"""
function make_weighted_raster_graph(
        weight_raster::GeoArray,
        vertex_raster::GeoArray;
        directed::Bool = false,
        condition_raster::GeoArray = weight_raster,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        connect_using_avg_weights::Bool = true,
        combine::Function = min
    )
    not_no_data = weight_raster .!= weight_raster.missingval
    if sum((weight_raster .<= 0) .& not_no_data) != 0
        @error("weight_raster contains 0 or negative values, " * 
               "which is not supported")
    end

    g = make_raster_graph(
        weight_raster,
        vertex_raster,
        directed = directed,
        weighted = true,
        condition_raster = condition_raster,
        condition = condition,
        cardinal_neighbors_only = cardinal_neighbors_only,
        connect_using_avg_weights = connect_using_avg_weights,
        combine = combine
    )

    return g
end


"""
    make_simple_raster_graph(
        raster::GeoArray,
        vertex_raster::GeoArray;
        directed::Bool = false,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        combine::Function = min
    )

Construct a `SimpleGraph` or `SimpleDiGraph` (if `directed = true`) based on 
two rasters (`raster` and `vertex_raster`). This function is useful 
if you already have a custom vertex raster and don't want SpatialGraphs.jl to
make one for you. The vertex raster denotes the spatial locations of each vertex
in the graph, and `raster` is used to construct the graph and determine which
vertices to connect.

## Parameters
`raster`: A GeoData.GeoArray on which to base the graph. Any pixel in `raster` 
with a value not equal to `raster.missingval` will be assigned a vertex
in the graph (corresponding to its centroid). The values in the raster can also
be used to determine which vertices to connect. See `condition` below for more
information.

`vertex_raster`: A `GeoData.GeoArray` with integer values ranging from 1:n, 
where n is the number of unique vertices in the graph. 

## Arguments
`directed`: A `Bool` determining whether the graph should be directed.

`condition`: A function that compares the values in `condition_raster` for two
neighbors to determine if those neighbors should be connected. The function must
compare two values and return either `true` or `false`. Useful functions to use 
here include `<`, `<=`, `==`, etc. The first argument to `condition` corresponds
to the source vertex, and the second argument corresponds to the destination 
vertex. So, if you only want to connect sources to destinations with a lower
value in `condition_raster` (e.g. in the case of developing a hydrologic flow 
graph based on elevation), then you would use `<` for `condition`. Defaults to 
`is_data`, which results in neighbors being connected as long as they are not 
NoData (`condition_raster.missingval`) in `condition_raster`. Note that if using
an inequality function (or any function where the result depends on argument 
position), the graph must be directed. For undirected graphs, you can use either
`is_data` or `==`, or any other custom function where argument position doesn't
matter, e.g. a function that determines whether the values in `condition_raster`
are within a certain distance of each other.

`cardinal_neighbors_only`: A `Bool` stating whether only cardinal neighbors
should be connected. By default, both cardinal _and_ diagonal neighbors are
connected. Note that when determining weights between diagonal neighbors, the
increased distance between them (as compared to the distance between cardinal
neighbors) is accounted for.

`connect_using_avg_raster_val`: `Bool`. This is intended to offer methods that
complement those used in Circuitscape.jl and Omniscape.jl. In this context,
weights (the values in `weight_raster`) are in units of electrical resistance. 
If `false`, the weight between two nodes with resistances R1 and R2 is 
calculated by converting resistance to conductances, taking the average, then 
taking the inverse of the result to convert back to resistance: 
`1 / ((1/R1 + 1/R2) / 2)`. `connect_using_avg_weights = false` correspondes to 
the default settings in Circuitscape. Defaults to `true', in which case the 
simple average of the weights (adjusted for distance in the case of diagonal 
neighbors) in `weight_raster` are used.

`combine`: In the case that there are multiple edges defined for a single pair
of vertices, how should the weight be chosen? Defaults to `min`. See the docs
for `SparseArrays.sparse()` for more information.
"""
function make_simple_raster_graph(
        raster::GeoArray,
        vertex_raster::GeoArray;
        directed::Bool = false,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        combine::Function = min
    )
    g = make_raster_graph(
        raster,
        vertex_raster,
        directed = directed,
        weighted = false,
        condition_raster = raster,
        condition = condition,
        cardinal_neighbors_only = cardinal_neighbors_only,
        combine = combine
    )

    return g
end


function make_raster_graph(
        raster::GeoArray,
        vertex_raster::GeoArray;
        directed::Bool = false,
        weighted::Bool = true,
        condition_raster::GeoArray = raster,
        condition::Function = is_data,
        cardinal_neighbors_only::Bool = false,
        connect_using_avg_weights::Bool = true,
        combine::Function = min
    )
    # Which averaging function to use
    card_avg = connect_using_avg_weights ? res_cardinal_avg : cond_cardinal_avg
    diag_avg = connect_using_avg_weights ? res_diagonal_avg : cond_diagonal_avg
    dims = size(raster)
    no_data_val = raster.missingval

    sources = Vector{Int64}()
    destinations = Vector{Int64}()
    weighted && (node_weights = Vector{Float64}())
    
    # Add the edges
    # Only need to do neighbors down or to the right for undirected graphs
    # because edge additions will be redundant.
    # Cardinal directions are in quotes since sometimes GeoArray are permuted
    # such that the top row doesn't necessarily correspond to the northern-most
    # pixels
    for row in 1:dims[1]
        for column in 1:dims[2]
            source_idx = CartesianIndex((row, column))
            if vertex_raster[source_idx] == 0 || 
                    raster[source_idx] == no_data_val
                continue
            else
                # Cardinal destination indices needed for undirected graph
                east_idx = CartesianIndex((row, column + 1))
                south_idx = CartesianIndex((row + 1, column))

                ## Add cardinal neighbors
                # "East"
                if column != dims[2] && vertex_raster[east_idx] != 0 && 
                        raster[east_idx] != no_data_val
                    if condition(condition_raster[source_idx],
                                 condition_raster[east_idx])
                        push!(sources, vertex_raster[source_idx])
                        push!(destinations, vertex_raster[east_idx])
                        if weighted
                            res = card_avg(raster[source_idx],
                                           raster[east_idx])
                            push!(node_weights, res)
                        end
                    end
                end
                # "South"
                if row != dims[1] && vertex_raster[south_idx] != 0  && 
                        raster[south_idx] != no_data_val
                    if condition(condition_raster[source_idx],
                                 condition_raster[south_idx])
                        push!(sources, vertex_raster[source_idx])
                        push!(destinations, vertex_raster[south_idx])
                        if weighted
                            res = card_avg(raster[source_idx],
                                       raster[south_idx])
                            push!(node_weights, res)
                        end
                    end
                end

                ## Add diagonal neighbors if needed
                if !cardinal_neighbors_only
                    # Diagonal destination indices needed for undirected graph
                    ne_idx = CartesianIndex((row - 1, column + 1))
                    se_idx = CartesianIndex((row + 1, column + 1))

                    # "Northeast"
                    if column != dims[2] && row != 1 && 
                            vertex_raster[ne_idx] != 0 && 
                            raster[ne_idx] != no_data_val
                        if condition(condition_raster[source_idx],
                                     condition_raster[ne_idx])
                            push!(sources, vertex_raster[source_idx])
                            push!(destinations, vertex_raster[ne_idx])
                            if weighted
                                res = diag_avg(raster[source_idx],
                                               raster[ne_idx])
                                push!(node_weights, res)
                            end
                        end
                    end
                    # "Southeast"
                    if row != dims[1] && column != dims[2] && 
                            vertex_raster[se_idx] != 0  && 
                            raster[se_idx] != no_data_val
                        if condition(condition_raster[source_idx],
                                     condition_raster[se_idx])
                            push!(sources, vertex_raster[source_idx])
                            push!(destinations, vertex_raster[se_idx])
                            if weighted
                                res = diag_avg(raster[source_idx],
                                           raster[se_idx])
                                push!(node_weights, res)
                            end
                        end
                    end
                end

                ## If it's a directed graph we also need to account for west,
                ## north, southwest, and northwest destinations
                if directed
                    # Cardinal neighbor destination indices for DiGraph
                    west_idx = CartesianIndex((row, column - 1))
                    north_idx = CartesianIndex((row - 1, column))

                    # "West"
                    if column != 1 && vertex_raster[west_idx] != 0 && 
                        raster[west_idx] != no_data_val
                        if condition(condition_raster[source_idx],
                                     condition_raster[west_idx])
                            push!(sources, vertex_raster[source_idx])
                            push!(destinations, vertex_raster[west_idx])
                            if weighted
                                res = card_avg(raster[source_idx],
                                           raster[west_idx])
                                push!(node_weights, res)
                            end
                        end
                    end

                    # "North"
                    if row != 1 && vertex_raster[north_idx] != 0 && 
                        raster[north_idx] != no_data_val
                        if condition(condition_raster[source_idx],
                                     condition_raster[north_idx])
                            push!(sources, vertex_raster[source_idx])
                            push!(destinations, vertex_raster[north_idx])
                            if weighted
                                res = card_avg(raster[source_idx],
                                               raster[north_idx])
                                push!(node_weights, res)
                            end
                        end
                    end

                    ## Diagonal neighbors for directed graphs
                    if !cardinal_neighbors_only
                        # Diagonal neighbor destination indices for DiGraph
                        nw_idx = CartesianIndex((row - 1, column - 1))
                        sw_idx = CartesianIndex((row + 1, column - 1))

                        # "Northwest"
                        if column != 1 && row != 1 && 
                            vertex_raster[nw_idx] != 0 && 
                            raster[nw_idx] != no_data_val
                            if condition(condition_raster[source_idx],
                                        condition_raster[nw_idx])
                                push!(sources, vertex_raster[source_idx])
                                push!(destinations, vertex_raster[nw_idx])
                                if weighted
                                    res = diag_avg(raster[source_idx],
                                                   raster[nw_idx])
                                    push!(node_weights, res)
                                end
                            end
                        end

                        # "Southwest"
                        if row != dims[1] && column != 1 && 
                            vertex_raster[sw_idx] != 0  && 
                            raster[sw_idx] != no_data_val
                            if condition(condition_raster[source_idx],
                                        condition_raster[sw_idx])
                                push!(sources, vertex_raster[source_idx])
                                push!(destinations, vertex_raster[sw_idx])
                                if weighted
                                    res = diag_avg(raster[source_idx],
                                               raster[sw_idx])
                                    push!(node_weights, res)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # push!'ing to vectors then combining is way faster than add_edge!
    # Sometimes, depending on the nodemap, there may be multiple edges with
    # different weights connecting the same two vertices (only when vertices
    # span multiple "pixels"). Default to using the min of all duplicates via
    # combine = min
    if directed
        if weighted
            g = SimpleWeightedDiGraph(
                sources,
                destinations,
                node_weights,
                combine = combine
            )
        else
            n_nodes = max(maximum(sources), maximum(destinations))
            g = SimpleDiGraph(
                sparse(sources, destinations, 1, n_nodes, n_nodes, combine)
            )
        end
    else
        if weighted
            g = SimpleWeightedGraph(
                sources,
                destinations,
                node_weights,
                combine = combine
            )
        else
            n_nodes = max(maximum(sources), maximum(destinations))
            g = SimpleGraph(
                sparse(
                    vcat(sources, destinations),
                    vcat(destinations, sources),
                    1,
                    n_nodes,
                    n_nodes,
                    combine
                )
            )
        end

    end

    return g
end

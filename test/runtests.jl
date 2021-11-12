using Rasters, Graphs, SimpleWeightedGraphs, SpatialGraphs, Test

@testset "Simple Raster Graph Construction" begin
    include("rastergraphs.jl")
end

@testset "Weighted Raster Graph Construction" begin
    include("weightedrastergraphs.jl")
end

@testset "Simple Raster DiGraph Construction" begin
    include("rasterdigraphs.jl")
end

@testset "Weighted Raster DiGraph Construction" begin
    include("weightedrasterdigraphs.jl")
end

@testset "Graphs Interface" begin
    include("graph_interface.jl")
end

printstyled("Checking that Base.show works...\n", bold = true)

A_array = Array{Float64}(undef, (3, 4, 1))
A_array[:,:,:] = [1, 3, 2, 0.5, 10, 8, 5, -9999, 3, 1, 2, 6]
x = X(1:4)
y = Y(1:3)
band = Band(1:1)

weight_raster = Raster(A_array, (y, x, band), missingval = -9999)
rasgraph = weightedrastergraph(weight_raster)
show(rasgraph);print("\n")


using Test, SpatialGraphs

@testset "Raster Graph Construction" begin
    include("rastergraphs.jl")
end

@testset "Raster DiGraph Construction" begin
    include("rasterdigraphs.jl")
end
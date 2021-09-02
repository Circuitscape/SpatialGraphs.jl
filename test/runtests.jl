using Test, SpatialGraphs

@testset "Weighted Raster Graph Construction" begin
    include("weightedrastergraphs.jl")
end

@testset "Simple Raster DiGraph Construction" begin
    include("simplerasterdigraphs.jl")
end

@testset "Weighted Raster DiGraph Construction" begin
    include("weightedrasterdigraphs.jl")
end
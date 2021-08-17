# Fixes GR warnings in Travis
ENV["GKSwstype"] = "100"

using Documenter, SpatialGraphs

const formats = Any[
    Documenter.HTML(
    	assets = [
    		"assets/custom.css"
    	],
	edit_link = :commit,
    ),
]

makedocs(
	format = formats,
    modules = [SpatialGraphs],
    authors = "Vincent A. Landau",
    sitename = "SpatialGraphs.jl",
    pages = ["About" => "index.md"],
)

deploydocs(
    repo = "github.com/Circuitscape/SpatialGraphs.jl.git",
    devbranch = "main",
    devurl = "latest"
)
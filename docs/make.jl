using Documenter, SpatialGraphs

const formats = Any[
    Documenter.HTML(
	    edit_link = :commit,
    )
]

makedocs(
	format = formats,
    modules = [SpatialGraphs],
    authors = "Vincent A. Landau",
    sitename = "SpatialGraphs.jl",
    pages = [
        "About" => "index.md",
        "Graph Types" => "graphtypes.md",
        "User Guide" => "userguide.md",
        "Examples" => "examples.md"
    ]
)

deploydocs(
    repo = "github.com/Circuitscape/SpatialGraphs.jl.git",
    devbranch = "main",
    devurl = "latest"
)
import Base: show

function show(io::IO, g::AbstractRasterGraph)
    printstyled("$(typeof(g))\n", color=:blue)
    printstyled("  graph", color=:red)
    print(": $(typeof(g.graph)) with $(nv(g)) vertices and $(ne(g)) edges\n")

    printstyled("  vertex_raster", color=:red)
    print(
        ": $(string(nameof(typeof(g.vertex_raster))))" *
        "{$(eltype(g.vertex_raster)), $(length(dims(g.vertex_raster)))}"
    )
    printstyled(" with dimensions:\n", color=:light_black)

    x_dim = dims(g.vertex_raster, X)
    y_dim = dims(g.vertex_raster, Y)
    printstyled("    X", color=:cyan)
    print(
        ": range($(minimum(x_dim)), $(maximum(x_dim)), step=$(x_dim.val[2] - x_dim.val[1]))\n"

    )
    printstyled("    Y", color=:cyan)
    print(
        ": range($(minimum(y_dim)), $(maximum(y_dim)), step=$(y_dim.val[2] - y_dim.val[1]))"
    )
    if (hasdim(g.vertex_raster, Band))
        band_dim = dims(g.vertex_raster, Band)
        printstyled("\n    Band", color=:cyan)
        print(
            ": $(minimum(band_dim)):$(maximum(band_dim))"
        )
    end
end
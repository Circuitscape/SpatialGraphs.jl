import Base: show

function show(io::IO, g::AbstractRasterGraph)
    printstyled(io, "$(typeof(g))", color=:blue)
    print(io, ":\n")
    printstyled(io, "  graph", color=:red)
    print(io, ": $(typeof(g.graph)) with $(nv(g)) vertices and $(ne(g)) edges\n")

    printstyled(io, "  vertex_raster", color=:red)
    print(
        io,
        ": $(string(nameof(typeof(g.vertex_raster))))" *
        "{$(eltype(g.vertex_raster)), $(length(dims(g.vertex_raster)))}"
    )
    printstyled(io, " with dimensions:\n", color=:light_black)

    x_dim = dims(g.vertex_raster, X)
    y_dim = dims(g.vertex_raster, Y)
    printstyled(io, "    X", color=:cyan)
    print(
        io,
        ": range($(minimum(x_dim)), $(maximum(x_dim)), step=$(x_dim.val[2] - x_dim.val[1]))\n"

    )
    printstyled(io, "    Y", color=:cyan)
    print(
        io,
        ": range($(minimum(y_dim)), $(maximum(y_dim)), step=$(y_dim.val[2] - y_dim.val[1]))"
    )
    if (hasdim(g.vertex_raster, Band))
        band_dim = dims(g.vertex_raster, Band)
        printstyled(io, "\n    Band", color=:cyan)
        print(
            io,
            ": $(minimum(band_dim)):$(maximum(band_dim))"
        )
    end
end
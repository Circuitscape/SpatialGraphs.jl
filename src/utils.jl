
"""
    is_data(a::Number, b::Number)

This function is the default used for deciding whether to connect two neighbors
in a GeoArray when constructing a graph. It always returns `true`, so all
neighbors will be connected

Returns `true`
"""
function is_data(a::Number, b::Number)
    return true
end

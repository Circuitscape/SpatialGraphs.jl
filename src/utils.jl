
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


# averaging functions, to operate on resistances
# cond_* functions convert to conductance, calulate mean, then convert to resistance
# res_* functions directly operate on the resistances
cond_cardinal_avg(x, y) = 1 / ((1/x + 1/y) / 2)
cond_diagonal_avg(x, y) = 1 / ((1/x + 1/y) / (2 * √2))
res_cardinal_avg(x, y) = (x + y) / 2
res_diagonal_avg(x, y) = ((x + y) * √2) / 2


# function get_dims(A::GeoData.GeoArray)
#     y_first = dims(mytif)[1] isa XDim
#     first_dim_type = y_first ? YDim : XDim
#     second_dim_type = y_first ? XDim : YDim

#     (dims(A, first_dim_type), dims(A, second_dim_type), Band(1:1, Categorical(order = Ordered())))
# end
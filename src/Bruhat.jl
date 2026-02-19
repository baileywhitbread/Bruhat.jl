module Bruhat

using Reexport
@reexport using Chevie
@reexport using Logging
@reexport using Random

include("intersections.jl")

export intersections_rational
export intersections_geometric

end

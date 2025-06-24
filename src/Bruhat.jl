module Bruhat

using Reexport
@reexport using Chevie
@reexport using Logging
@reexport using Random

export intersections_rational
export intersections_geometric

include("intersections.jl")

end

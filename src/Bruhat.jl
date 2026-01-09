module Bruhat

using Reexport
@reexport using Chevie
@reexport using Logging
@reexport using Random

include("checks.jl")
include("intersections.jl")

export intersections_rational
export intersections_geometric
export intersections_geometric_ordered
export verify_lusztig

end

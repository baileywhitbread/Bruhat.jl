module Bruhat

using Reexport
@reexport using Chevie
@reexport using Logging
@reexport using Random

export intersections

include("intersections.jl")

end

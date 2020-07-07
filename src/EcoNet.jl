module EcoNet

# dependency
using LinearAlgebra, Statistics, Distributions

export
# foodwebs.jl
FoodWeb, generateFoodWeb,

# odeModel.jl
webDynamic!

include("foodwebs.jl")
include("odeModel.jl")

end # module

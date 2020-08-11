module EcoNet

# dependency
using LinearAlgebra, Statistics, Distributions

export
# foodwebs.jl
FoodWeb, generateFoodWeb, info,

# odeModel.jl
webDynamic!

include("foodwebs.jl")
include("odeModel.jl")

end # module

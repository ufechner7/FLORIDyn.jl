module FLORIDyn

using Interpolations, LinearAlgebra

export Direction_Constant, Direction_Constant_wErrCov
export getWindDirT, getWindDirT_EnKF

# the different wind direction types
struct Direction_Constant end
struct Direction_Constant_wErrCov end

include("windfield.jl")

end

module FLORIDyn

using Interpolations, LinearAlgebra

export Direction_Constant, Direction_Constant_wErrCov, Direction_EnKF_InterpTurbine
export getWindDirT, getWindDirT_EnKF

# the different wind direction types (dir_mode)
struct Direction_Constant end
struct Direction_Constant_wErrCov end
struct Direction_EnKF_InterpTurbine end
struct Direction_Interpolation end

include("windfield.jl")

end

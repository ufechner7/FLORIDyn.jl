module FLORIDyn

using Interpolations

export Direction_Constant
export getWindDirT, getWindDirT_EnKF

# define your different modes for the function
struct Direction_Constant end
struct Mode2 end

include("windfield.jl")

end

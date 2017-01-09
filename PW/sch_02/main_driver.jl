# These must be set manually :-(
include("../common/PWGrid_v02.jl")
#include("../common/PWGrid_v03.jl")

include("driver_harm.jl")

@time test_main( Ns_in=[30, 30, 30] )
#@time test_main( ecutwfc_Ry=54.0 )

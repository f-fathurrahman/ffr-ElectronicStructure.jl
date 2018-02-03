include("pulaymix.jl")

function main_debug()
    Npoints = 100
    rho = rand(Npoints)
    rhonew = rand(Npoints)
    beta = 0.5
    MIXDIM = 4
    dvmat = rand(Float64,Npoints,MIXDIM)
    vmat  = rand(Float64,Npoints,MIXDIM)
    iter = 5
    rho = pulaymix!( rho, rhonew, beta, dvmat, vmat, iter, MIXDIM)
end

main_debug()
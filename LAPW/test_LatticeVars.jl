include("LatticeVars.jl")

function main()
    LatVecs = zeros(3,3)
    LatVecs[1,:] = [5.0, 5.0, 0.0]
    LatVecs[2,:] = [5.0, 0.0, 5.0]
    LatVecs[3,:] = [0.0, 6.0, 5.0]
    
    lattice_vars = LatticeVars( LatVecs )
    avec = lattice_vars.avec
    bvec = lattice_vars.bvec
    res = (avec*bvec')/(2*pi)
    display(res); println()

    println("Pass here")
end

main()
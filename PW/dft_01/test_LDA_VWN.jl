include("LDA_VWN.jl")

function test_main()
    rhoe = 1.5*ones(5)
    epsxc = excVWN( rhoe )
    depsxc = excpVWN( rhoe )
    println(epsxc)
    println(depsxc)
end

test_main()

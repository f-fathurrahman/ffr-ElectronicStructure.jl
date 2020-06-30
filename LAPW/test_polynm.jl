using Printf

include("polynm.jl")

function my_sin(x)
    return sin(10*x)
end

function integ_my_sin(x)
    return -cos(10*x)/10
end

function d_my_sin(x)
    return 10*cos(10*x)
end

function test_main()
    Npoints = 10
    xa = range(1.0, 1.5, length=Npoints)
    ya = my_sin.(xa)

    m = 0
    for i in 1:Npoints
        x = xa[i]
        y = polynm(m, Npoints, xa, ya, x)
        @printf("%18.10f %18.10f %18.10f\n", x, y, ya[i])
    end

    x_plt = range(1.0, 1.5, length=5*Npoints)
    for x in x_plt
        y = polynm(m, Npoints, xa, ya, x)
        #@printf("%18.10f %18.10f %18.10f\n", x, y, my_sin(x))
        #@printf("%18.10f %18.10f %18.10f\n", x, y, integ_my_sin(x))
    end
end

test_main()
struct PGBF
    center::Tuple3F64
    power::Tuple3I64
    expn::Float64
    NORM::Float64
end

import Base.println
function println( bf::PGBF )
    @printf("\n")
    @printf("Primitive Gaussian:\n")
    @printf("Exponent: %18.10f\n", bf.expn)
    @printf("Center  : (%18.10f,%18.10f,%18.10f)\n", bf.center[1], bf.center[2], bf.center[3])
    @printf("AM      : (%2d,%2d,%2d)\n", bf.power[1], bf.power[2], bf.power[3])
    @printf("Norm    : %18.10f\n", bf.NORM)
    @printf("\n")
end

function PGBF( expn::Float64, center::Tuple3F64, power::Tuple3I64 )
    NORM = 1.0 / sqrt(overlap(expn,center,power))
    p = PGBF( center, power, expn, NORM )
    return p
end

function PGBF( expn::Float64 )
    center = (0.0,0.0,0.0)
    power = (0,0,0)
    NORM = 1.0 / sqrt(overlap(expn,center,power))
    p = PGBF( center, power, expn, NORM )
    return p
end

function evaluate( bf::PGBF, x::Float64, y::Float64, z::Float64 )
    xo = bf.center[1]
    yo = bf.center[2]
    zo = bf.center[3]
    I = bf.power[1]
    J = bf.power[2]
    K = bf.power[3]
    expn = bf.expn
    NORM = bf.NORM

    dx, dy, dz = x-xo, y-yo, z-zo
    r2 = dist2(dx, dy, dz)
    return NORM*(dx^I)*(dy^J)*(dz^K)*exp(-expn*r2)
end


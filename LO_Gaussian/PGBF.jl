struct PGBF
    expn::Float64
    center::Tuple3F64
    power::Tuple3I64
    norm::Float64
end

import Base.println
function println( bf::PGBF )
    @printf("\n")
    @printf("Info for PGBF\n")
    @printf("Exponent: %18.10f\n", bf.expn)
    @printf("Center: (%18.10f,%18.10f,%18.10f)\n", bf.center[1], bf.center[2], bf.center[3])
    @printf("Angular momentum: (%2d,%2d,%2d)\n", bf.power[1], bf.power[2], bf.power[3])
    @printf("Norm: %18.10f\n", bf.norm)
    @printf("\n")
end

function PGBF( expn, center::Tuple3F64, power::Tuple3I64 )
    norm = 1.0 / sqrt(overlap(expn,center,power))
    p = PGBF( expn, center, power, norm )
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

    dx, dy, dz = x-xo, y-yo, z-zo
    r2 = dist2(dx, dy, dz)
    return bf.norm*(dx^I)*(dy^J)*(dz^K)*exp(-expn*r2)
end


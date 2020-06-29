# TODO: make this immutable, possibly by changing the mechanism to build
# the `BasisSet`

"""
Contracted Gaussian basis function.
```
CGBF = sum_{i=Ncontr} coefs[i]*PGBF(expn, center, power)
```
"""
mutable struct CGBF
    center::Tuple3F64
    power::Tuple3I64
    pgbfs::Array{PGBF,1}
    coefs::Array{Float64,1}
    NORM::Float64
end

import Base.println
function println( bf::CGBF )
    Ncontr = length(bf.pgbfs)
    @printf("\nContracted Gaussian basis function\n")
    @printf("\nNcontr = %d\n\n", Ncontr)
    for i = 1:Ncontr
        @printf("Coef = %18.10f", bf.coefs[i])
        println(bf.pgbfs[i])
    end
end

function init_CGBF( center::Tuple3F64, power::Tuple3I64 )
    return CGBF( center, power, PGBF[], Float64[], 1.0 )
end

init_CGBF(x=0,y=0,z=0,I=0,J=0,K=0) = CGBF( (x,y,z), (I,J,K), PGBF[], Float64[], 1.0 )

function evaluate(bf::CGBF,x::Float64,y::Float64,z::Float64)
    s = 0
    for (c,pbf) in primitives(bf)
        s += c*evaluate(pbf,x,y,z)
    end
    return bf.NORM*s
end

function normalize!(bf::CGBF)
    bf.NORM /= sqrt(overlap(bf,bf))
end

primitives(a::CGBF) = zip(a.coefs,a.pgbfs)

# This should become the constructor along with center and power
import Base: push!
function push!(cbf::CGBF, expn, coef)
    Base.push!(cbf.pgbfs, PGBF(expn, cbf.center, cbf.power))
    Base.push!(cbf.coefs, coef)
    normalize!(cbf)
end

function contract(f, a::CGBF, b::CGBF)
    s = 0.0
    for (ca,abf) in primitives(a)
    for (cb,bbf) in primitives(b)
        s += ca*cb*f(abf,bbf)
    end
    end
    return a.NORM*b.NORM*s
end

function contract(f,a::CGBF,b::CGBF,c::CGBF,d::CGBF)
    s = 0
    for (ca,abf) in primitives(a)
    for (cb,bbf) in primitives(b)
    for (cc,cbf) in primitives(c)
    for (cd,dbf) in primitives(d)
        s += ca*cb*cc*cd*f(abf,bbf,cbf,dbf)
    end
    end
    end
    end
    return a.NORM*b.NORM*c.NORM*d.NORM*s
end

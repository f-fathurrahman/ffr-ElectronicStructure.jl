### Nuclear attraction term

function Aterm( i::Int64, r::Int64, u::Int64, l1::Int64, l2::Int64,
                ax::Float64,bx::Float64,cx::Float64, gamma::Float64)
    term1 = (-1)^(i+u)*binomial_prefactor(i,l1,l2,ax,bx)
    term2 = factorial(i)*cx^(i-2r-2u)
    term3 = (1/4/gamma)^(r+u)/factorial(r)/factorial(u)/factorial(i-2r-2u)
    return term1*term2*term3
end

function Aarray(l1::Int64,l2::Int64,a::Float64,b::Float64,c::Float64,g::Float64)
    Imax = l1 + l2 + 1
    A = zeros(Float64,Imax)
    for i in 0:(Imax-1)
        for r in 0:div(i,2)
            for u in 0:div(i-2r,2)
                I = i-2r-u+1
                A[I] += Aterm(i,r,u,l1,l2,a,b,c,g)
            end
        end
    end
    return A
end

function nuclear_attraction(
        aexpn::Float64, acenter::Tuple3F64, apower::Tuple3I64,
        bexpn::Float64, bcenter::Tuple3F64, bpower::Tuple3I64,
        ccenter::Tuple3F64
    )
    
    ax, ay, az = acenter[1],acenter[2],acenter[3]
    bx, by, bz = bcenter[1],bcenter[2],bcenter[3]
                  
    aI,aJ,aK = apower[1],apower[2],apower[3]
    bI,bJ,bK = bpower[1],bpower[2],bpower[3]

    cx, cy, cz = ccenter[1],ccenter[2],ccenter[3]

    px,py,pz = gaussian_product_center(aexpn,ax,ay,az,bexpn,bx,by,bz)
    gamma = aexpn + bexpn
    
    rab2 = dist2(ax-bx,ay-by,az-bz)
    rcp2 = dist2(cx-px,cy-py,cz-pz)
    
    Ax = Aarray(aI,bI,px-ax,px-bx,px-cx,gamma)
    Ay = Aarray(aJ,bJ,py-ay,py-by,py-cy,gamma)
    Az = Aarray(aK,bK,pz-az,pz-bz,pz-cz,gamma)
    total = 0
    for I in 0:(aI+bI)
        for J in 0:(aJ+bJ)
            for K in 0:(aK+bK)
                total += Ax[I+1]*Ay[J+1]*Az[K+1]*Fgamma(I+J+K,rcp2*gamma)
            end
        end
    end
    val = -2.0*pi*exp(-aexpn*bexpn*rab2/gamma)*total/gamma
    return val
end

function nuclear_attraction( a::PGBF, b::PGBF, ccenter::Tuple3F64, atno )
    Vab = nuclear_attraction(a.expn, a.center, a.power,
                             b.expn, b.center, b.power, ccenter)
    return atno*a.NORM*b.NORM*Vab
end

#nuclear_attraction(a::PGBF,b::PGBF,c::Atom) = c.atno*nuclear_attraction(a,b,c.x,c.y,c.z)

function nuclear_attraction( a::PGBF, b::PGBF, atoms::Atoms )
    Natoms = atoms.Natoms
    Vab = 0.0
    Zatoms = get_Zatoms( atoms )
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        cx = atoms.positions[1,ia]
        cy = atoms.positions[2,ia]
        cz = atoms.positions[3,ia]
        isp = atm2species[ia]
        atno = Zatoms[isp]
        Vab = Vab + nuclear_attraction( a,b,(cx,cy,cz),atno )
    end
    return Vab
end

function Fgamma(m::Int64,x::Float64,SMALL::Float64=1e-12)
    x = max(x,SMALL) # Evidently needs underflow protection
    return 0.5*x^(-m-0.5)*gammainc(m+0.5,x)
end

function gammainc(a::Float64,x::Float64)
    # This is the series version of gamma from pyquante. For reasons I don't get, it
    # doesn't work around a=1. This works alright, but is only a stopgap solution
    # until Julia gets an incomplete gamma function programmed
    if abs(a-1) < 1e-3
        println("Warning: gammainc_series is known to have problems for a ~ 1")
    end
    if x < (a+1.0)
        #Use the series representation
        gam,gln = gser(a,x)
    else
        #Use continued fractions
        gamc,gln = gcf(a,x)
        gam = 1-gamc
    end
    return exp(gln)*gam
end

function gser( a::Float64, x::Float64, ITMAX::Int64=100, EPS::Float64=3e-9 )
    # Series representation of Gamma. NumRec sect 6.1.
    gln = lgamma(a)
    if x == 0
        return 0,gln
    end
    ap = a
    delt = s = 1/a
    for i in 1:ITMAX
        ap += 1
        delt *= (x/ap)
        s += delt
        if abs(delt) < abs(s)*EPS
            break
        end
    end
    return s*exp(-x+a*log(x)-gln),gln
end

function gcf(a::Float64,x::Float64,ITMAX::Int64=200,EPS::Float64=3e-9,FPMIN::Float64=1e-30)
    # Continued fraction representation of Gamma. NumRec sect 6.1"
    gln = lgamma(a)
    b = x + 1.0 - a
    c = 1.0/FPMIN
    d = 1.0/b
    h = d
    for i in 1:ITMAX
        an = -i*(i-a)
        b = b + 2.0
        d = an*d + b
        if abs(d) < FPMIN
            d = FPMIN
        end
        c = b + an/c
        if abs(c) < FPMIN
            c = FPMIN
        end
        d = 1.0/d
        delt = d*c
        h = h*delt
        if abs(delt-1.0) < EPS
            break
        end
    end
    gammcf = exp(-x + a*log(x) - gln)*h
    return gammcf,gln
end


# Need a nested scope to squeeze this into the contract function
function nuclear_attraction( a::CGBF, b::CGBF, ccenter::Tuple3F64, atno )
    na(a,b) = nuclear_attraction( a, b, ccenter, atno )
    contract(na,a,b)
end

#function nuclear_attraction(a::CGBF,b::CGBF,c::Atom)
#  na(a,b) = nuclear_attraction(a,b,c)
#  contract(na,a,b)
#end

function nuclear_attraction( a::CGBF, b::CGBF, atoms::Atoms)
    na(a,b) = nuclear_attraction( a, b, atoms )
    contract(na,a,b)
end


function nuclear_repulsion( ZA, posA, ZB, posB )
    a.atno*b.atno/sqrt( dist2( ax-bx, ay-by, az-bz) )
end

function nuclear_repulsion( atoms::Atoms )
    nr = 0.0
    Zatoms = get_Zatoms( atoms )
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        for ja = ia+1:Natoms
            isp = atm2species[ia]
            jsp = atm2species[ja]
            Zi = Zatoms[isp]
            Zj = Zatoms[jsp]
            dr = atoms.positions[:,ia] - atoms.positions[:,ja]
            nr = nr + Zi*Zj/norm(dr)
        end
    end
    return nr
end

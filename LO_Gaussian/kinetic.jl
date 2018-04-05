"""
### Kinetic matrix elements
"""
function kinetic( a::PGBF, b::PGBF )
    K = kinetic(a.expn, a.center, a.power, b.expn, b.center, b.power)
    return a.NORM*b.NORM*K
end

function kinetic( aexpn::Float64, acenter::Tuple3F64, apower::Tuple3I64,
                  bexpn::Float64, bcenter::Tuple3F64, bpower::Tuple3I64 )
    
    ax, ay, az = acenter[1],acenter[2],acenter[3]
    bx, by, bz = bcenter[1],bcenter[2],bcenter[3]
              
    aI,aJ,aK = apower[1],apower[2],apower[3]
    bI,bJ,bK = bpower[1],bpower[2],bpower[3]

    overlap0 = overlap( aexpn, acenter, apower, bexpn, bcenter, bpower )
    
    overlapx1 = overlap( aexpn, acenter, apower, bexpn, bcenter, (bI+2,bJ,bK) )
    
    overlapy1 = overlap( aexpn, acenter, apower, bexpn, bcenter, (bI,bJ+2,bK) )
    
    overlapz1 = overlap( aexpn, acenter, apower, bexpn, bcenter, (bI,bJ,bK+2) )
    
    overlapx2 = overlap( aexpn, acenter, apower, bexpn, bcenter, (bI-2,bJ,bK) )
    
    overlapy2 = overlap( aexpn, acenter, apower, bexpn, bcenter, (bI,bJ-2,bK) )
    
    overlapz2 = overlap( aexpn, acenter, apower, bexpn, bcenter, (bI,bJ,bK-2) )
  
    term0 = bexpn*( 2 * (bI + bJ + bK) + 3)*overlap0
    term1 = -2*(bexpn^2)*(overlapx1 + overlapy1 + overlapz1)
    term2 = -0.5*(bI*(bI-1)*overlapx2+bJ*(bJ-1)*overlapy2 + bK*(bK-1)*overlapz2)
    return term0 + term1 + term2
end

kinetic( a::CGBF, b::CGBF ) = contract( kinetic, a, b )

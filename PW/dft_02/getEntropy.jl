function getEntropy(occ,Tbeta)
    const SMALL = 1.e-10
    # nonzero occupations
    inz = find(occ > SMALL)
    occnz = occ(inz)
    e = sum( occnz .* log(occnz) )
    #
    # zero occupations
    #
    iz = find(abs(occ) <= SMALL)
    occz = occ(iz)
    e = e + sum( (1.0 - occz) .* log( 1.0 - occz ) )
    e = e / Tbeta
    return e
end
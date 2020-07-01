function mtdmin(
    latt_vars::LatticeVars, atm_vars::AtomicVars, rmt
)
    v1 = zeros(Float64,3)
    v2 = zeros(Float64,3)
    is = 1
    js = 1
    dmin = 1.0e6

    avec = latt_vars.avec
    epslat = latt_vars.epslat
    atposc = atm_vars.atposc
    nspecies = atm_vars.nspecies
    natoms = atm_vars.natoms
    
    for i1 in -1:1, i2 in -1:1, i3 in -1:1
        #
        @views v1[:] = i1*avec[:,1] + i2*avec[:,2] + i3*avec[:,3]
        #
        for ks in 1:nspecies, ka in 1:natoms[ks]
            #
            @views v2[:] = v1[:] + atposc[:,ka,ks]
            #
            for ls in 1:nspecies
                #
                t1 = rmt[ks] + rmt[ls]
                #
                for la in 1:natoms[ls]
                    if ( (i1 != 0) || (i2 != 0) || (i3 != 0) || (ks != ls) ||  (ka != la) )
                        t2 = sqrt( (v2[1] - atposc[1,la,ls])^2 +
                                   (v2[2] - atposc[2,la,ls])^2 +
                                   (v2[3] - atposc[3,la,ls])^2 )
                        t3 = t2 - t1
                        if t3 < (dmin - epslat)
                            is = ks
                            js = ls
                            dmin = t3
                        end
                    end
                end
            end
        end
    end

    return dmin, is, js

end


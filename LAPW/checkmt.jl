function checkmt!(
    latt_vars::LatticeVars, atm_vars::AtomicVars, spsymb,
    mtr_vars::MuffinTinRadialVars
)

    rmt = mtr_vars.rmt
    rmtdelta = mtr_vars.rmtdelta
    epslat = latt_vars.epslat
    nspecies = atm_vars.nspecies

    rmt0 = zeros(Float64,nspecies)
    rmt0[1:nspecies] = rmt[1:nspecies]

    println("rmtdelta = ", rmtdelta)

    while true
        # find the minimum distance between muffin-tin surfaces
        dmin, is, js = mtdmin(latt_vars, atm_vars, rmt)
        # adjust muffin-tin radii if required1
        if dmin < (rmtdelta - epslat)
            println("Adjusting rmt")
            t1 = rmt[is] + rmt[js]
            t2 = (t1 + dmin - rmtdelta)/t1
            rmt[is] = rmt[is]*t2
            if is != js
                rmt[js] = rmt[js]*t2
            end
        else
            break
        end
    end

    for is in 1:nspecies
        if rmt[is] < 0.25
            @printf("Error(checkmt): muffin-tin radius too small for species %d %s\n", is, spsymb[is])
            @printf("Radius : %18.10f\n", rmt[is])
        end    
        if rmt[is] < rmt0[is]
            @printf("Info(checkmt): reduced muffin-tin radius of species %3d %s", is, spsymb[is])
            @printf(" is reduced from %f to %f\n", rmt0[is], rmt[is])
        end
    end

    return
end

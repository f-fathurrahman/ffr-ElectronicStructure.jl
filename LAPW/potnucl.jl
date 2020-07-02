function potnucl!(ptnucl::Bool, nr, r, zn, vn )
    SMALL = 1e-10
    if abs(zn) <= SMALL
        vn[:] .= 0.0
        return
    end
    if ptnucl
        # nucleus is taken to be a point particle
        for ir in 1:nr
            vn[ir] = zn/r[ir]
        end
    else
        # approximate nuclear radius
        rn = radnucl(zn)
        t1 = zn/(2.0*rn^3)
        t2 = 3.0*rn^2
        for ir in 1:nr
          if r[ir] < rn
            vn[ir] = t1*(t2 - r[ir]^2)
          else
            vn[ir] = zn/r[ir]
          end
        end
    end
    return
end

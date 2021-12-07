# !INPUT/OUTPUT PARAMETERS:
#   sol    : speed of light in atomic units (in,real)
#   ptnucl : .true. if the nucleus is a point particle (in,logical)
#   zn     : nuclear charge (in,real)
#   nst    : number of states to solve for (in,integer)
#   n      : priciple quantum number of each state (in,integer(nst))
#   l      : quantum number l of each state (in,integer(nst))
#   k      : quantum number k (l or l+1) of each state (in,integer(nst))
#   occ    : occupancy of each state (inout,real(nst))
#   xctype : exchange-correlation type (in,integer(3))
#   xcgrad : 1 for GGA functional, 0 otherwise (in,integer)
#   nr     : number of radial mesh points (in,integer)
#   r      : radial mesh (in,real(nr))
#   eval   : eigenvalue without rest-mass energy for each state (out,real(nst))
#   rho    : charge density (out,real(nr))
#   vr     : self-constistent potential (out,real(nr))
#   rwf    : major and minor components of radial wavefunctions for each state
#            (out,real(nr,2,nst))
function solve_atom!(
    sol, ptnucl, zn, nst, n, l, k, occ,
    xctype, xcgrad, nr, r, evals, rho, vr, rwf
)

# real(8), intent(in) :: sol
# logical, intent(in) :: ptnucl
# real(8), intent(in) :: zn
# integer, intent(in) :: nst
# integer, intent(in) :: n(nst),l(nst),k(nst)
# real(8), intent(inout) :: occ(nst)
# integer, intent(in) :: xctype(3),xcgrad
# integer, intent(in) :: nr
# real(8), intent(in) :: r(nr)
# real(8), intent(out) :: eval(nst)
# real(8), intent(out) :: rho(nr),vr(nr)
# real(8), intent(out) :: rwf(nr,2,nst)
  
    #local variables
    maxscl = 200

    # potential convergence tolerance
    SMALL = 1.0e-6

    @assert nst > 0

    # allocate local arrays
    vn = zeros(Float64,nr)
    vh = zeros(Float64,nr)
    ex = zeros(Float64,nr)
    ec = zeros(Float64,nr)
    vx = zeros(Float64,nr)
    vc = zeros(Float64,nr)
    vrp = zeros(Float64,nr)
    
    ri = zeros(Float64,nr)
    wpr = zeros(Float64,4,nr)
    fr1 = zeros(Float64,nr)
    fr2 = zeros(Float64,nr)
    gr1 = zeros(Float64,nr)
    gr2 = zeros(Float64,nr)
    
    #if (xcgrad.eq.1) then
    #  allocate(grho(nr),g2rho(nr),g3rho(nr))
    #end if

    # find total electronic charge
    ze = 0.0
    for ist in 1:nst
        ze = ze + occ[ist]
    end

    # set up nuclear potential
    potnucl!(ptnucl, nr, r, zn, vn)

    for ir in 1:nr
        ri[ir] = 1.0/r[ir]
        # initialise the Kohn-Sham potential to the nuclear potential
        vr[ir] = vn[ir]
    end

    # determine the weights for radial integration
    wsplintp!(nr, r, wpr)
    
    # initialise mixing parameter
    betamix = 0.5
    
    # initialise eigenvalues to relativistic values (minus the rest mass energy)
    for ist in 1:nst
        t1 = sqrt( k[ist]^2 - (zn/sol)^2 )
        t1 = ( n[ist] - abs(k[ist]) + t1)^2
        t1 = 1.0 + (zn/sol)^2/t1
        evals[ist] = sol^2/sqrt(t1) - sol^2
    end
    
    dvp = 0.0
    # start of self-consistent loop
    is_converged = false
    for iscl in 1:maxscl

        # solve the Dirac equation for each state
        for ist in 1:nst
            @views evals[ist] = rdirac!( sol, n[ist], l[ist], k[ist], nr, r, vr, evals[ist],
                rwf[:,1,ist], rwf[:,2,ist] )
        end
    
        # compute the charge density
        for ir in 1:nr
            ss = 0.0
            for ist in 1:nst
              ss = ss + occ[ist]*( rwf[ir,1,ist]^2 + rwf[ir,2,ist]^2 )
            end
            fr1[ir] = ss
            fr2[ir] = ss*ri[ir]
            rho[ir] = ss*ri[ir]^2 / (4*pi)
        end
        splintwp!(nr, wpr, fr1, gr1)
        splintwp!(nr, wpr, fr2, gr2)
        
        # find the Hartree potential
        t1 = gr2[nr]
        for ir in 1:nr
           vh[ir] = gr1[ir]*ri[ir] + t1 - gr2[ir]
        end
    
        # normalise charge density and potential
        t1 = ze/gr1[nr]
        rho[:] .= t1*rho[:]
        vh[:] .= t1*vh[:]
    
        # compute the exchange-correlation energy and potential
        if xcgrad == 1
            # GGA functional
            # |grad rho|
            fderiv!(1, nr, r, rho, grho)
            # grad^2 rho
            fderiv!(2, nr, r, rho, g2rho)
            for ir in 1:nr
                g2rho[ir] = g2rho[ir] + 2.0*ri[ir]*grho[ir]
            end
            # approximate (grad rho).(grad |grad rho|)
            for ir in 1:nr
                g3rho[ir] = grho[ir]*g2rho[ir]
            end
            #xcifc(xctype,n=nr,rho=rho,grho=grho,g2rho=g2rho,g3rho=g3rho,ex=ex, ec=ec,vx=vx,vc=vc)
        else
            # LDA functional
            #call xcifc(xctype,n=nr,rho=rho,ex=ex,ec=ec,vx=vx,vc=vc)
            for ir in 1:nr
                if rho[ir] <= 10e-12
                    #println("small rho = ", rho[ir])
                    ex[ir] = vx[ir] = 0.0
                    ec[ir] = vc[ir] = 0.0
                else
                    ex[ir], vx[ir] = XC_x_slater( rho[ir] )
                    ec[ir], vc[ir] = XC_c_vwn( rho[ir] )
                end
            end
        end
        
        #println("sum(vx) = ", sum(vx))
        #println("sum(vc) = ", sum(vc))
        #exit()
    
        # self-consistent potential
        @views vr[1:nr] = vh[1:nr] + vx[1:nr] + vc[1:nr]
    
        # determine change in potential
        ss = 0.0
        for ir in 1:nr
          ss = ss + ( vr[ir] - vrp[ir] )^2
        end
        dv = sqrt(ss)/nr
    
        if iscl > 2
            # reduce beta if change in potential is diverging
            if dv > dvp
                betamix = betamix*0.8
            end
            betamix = max(betamix,0.01)
        end
    
        dvp = dv
    
        for ir in 1:nr
            # mix old and new potentials
            vr[ir] = (1.0 - betamix)*vrp[ir] + betamix*vr[ir]
            vrp[ir] = vr[ir]
            # add nuclear potential
            vr[ir] = vr[ir] + vn[ir]
        end

        @printf("iscl = %5d dv = %18.10f\n", iscl, dv)

        # check for convergence
        if ( (iscl > 2) && (dv < SMALL) )
            println(".... Converged ....")
            is_converged = true
            break
        end
    end

    if !is_converged
        println("Warning(atom): maximum iterations exceeded")
    end

    return

end 


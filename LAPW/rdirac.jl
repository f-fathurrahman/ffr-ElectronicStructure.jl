# !INPUT/OUTPUT PARAMETERS:
#   sol  : speed of light in atomic units (in,real)
#   n    : principal quantum number (in,integer)
#   l    : quantum number l (in,integer)
#   k    : quantum number k (l or l+1) (in,integer)
#   nr   : number of radial mesh points (in,integer)
#   r    : radial mesh (in,real(nr))
#   vr   : potential on radial mesh (in,real(nr))
#   evals : eigenvalue without rest-mass energy (inout,real)
#   g0   : major component of the radial wavefunction (out,real(nr))
#   f0   : minor component of the radial wavefunction (out,real(nr))
# !DESCRIPTION:
#   Finds the solution to the radial Dirac equation for a given potential $v(r)$
#   and quantum numbers $n$, $k$ and $l$. The method involves integrating the
#   equation using the predictor-corrector method and adjusting $E$ until the
#   number of nodes in the wavefunction equals $n-l-1$. The calling routine must
#   provide an initial estimate for the eigenvalue. Note that the arrays
#   {\tt g0} and {\tt f0} represent the radial functions multiplied by $r$.

function rdirac!(sol, n, l, k, nr, r, vr, evals, g0, f0)

    # ! arguments
    # real(8), intent(in) :: sol
    # integer, intent(in) :: n,l,k,nr
    # real(8), intent(in) :: r(nr),vr(nr)
    # real(8), intent(inout) :: evals
    # real(8), intent(out) :: g0(nr),f0(nr)

    @assert k > 0
    @assert nr >= 4

    # local variables
    maxit = 2000
    # energy convergence tolerance
    SMALL = 1.e-12
    # automatic arrays
    g1 = zeros(Float64,nr)
    f1 = zeros(Float64,nr)
    fr = zeros(Float64,nr)

    if k > n
        println("n = ", n)
        println("k = ", k)
        error("Error incompatible n and k")
    end

    if ( (k == n) && (l != k-1) )
        println("n = ", n)
        println("k = ", k)
        println("l = ", k)
        error("Incompatible n, k and l")
    end

    if k == l
        kpa = k
    elseif k == (l+1)
        kpa = -k
    else
        println("l = ", l)
        println("k = ", k)
        println("Error: incompatible l and k")
    end

    #println("n   = ", n)
    #println("l   = ", l)
    #println("k   = ", k)
    #println("kpa = ", kpa)

    de = 1.0
    nndp = 0
    it = 0
    #
    while true
        #
        it = it + 1
        if it > maxit
            break
        end 
        # integrate the Dirac equation
        #println("calling rdiracint")
        #println("kpa = ", kpa)
        #println("nr  = ", nr)
        nn, evals = rdiracint!(sol, kpa, evals, nr, r, vr, g0, g1, f0, f1)
        # check the number of nodes
        nnd = nn - (n-l-1)
        if nnd > 0
            evals = evals - de
        else
            evals = evals + de
        end
        #
        if it > 1
            if ( (nnd != 0) || ( nndp != 0) )
                if nnd*nndp <= 0
                    de = de*0.5
                else
                    de = de*1.1
                end
            end
        end
        #
        nndp = nnd
        if de < SMALL*(abs(evals) + 1.0)
            break
        end
    end

    if it > maxit
        println("Warning(rdirac): maximum iterations exceeded")
    end

    # find effective infinity and set wavefunction to zero after that point
    # major component
    irm = nr
    for ir in 2:nr
        if ( (g0[ir-1]*g0[ir] < 0.0) || (g1[ir-1]*g1[ir] < 0.0) )
            irm = ir
        end
    end
    g0[irm:nr] .= 0.0
    # minor component
    irm = nr
    for ir in 2:nr
        if ( (f0[ir-1]*f0[ir] < 0.0) || (f1[ir-1]*f1[ir] < 0.0) )
            irm = ir
        end
    end
    f0[irm:nr] .= 0.0
    # normalise
    for ir in 1:nr
      fr[ir] = g0[ir]^2 + f0[ir]^2
    end

    t1 = splint(nr, r, fr)
    t1 = sqrt(abs(t1))

    if t1 > 0.0
      t1 = 1.0/t1
    else
      error("Error(rdirac): zero wavefunction")
    end
    g0[:] = t1*g0[:]
    f0[:] = t1*f0[:]
    return evals
end

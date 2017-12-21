function linsolve_pcg( Lmat::SparseMatrixCSC{Float64,Int64},
                       prec::SparseMatrixCSC{Float64,Int64},
                       b::Array{Float64,1};
                       x0 = nothing,
                       NiterMax = 1000, TOL=5.e-10,
                       convmsg=false, showprogress=false )
    #
    Npoints = size(b)[1]
    if x0 == nothing
        x = zeros( Float64, Npoints )
    else
        x = copy(x0)
    end
    #
    r = zeros( Float64, Npoints )
    p = zeros( Float64, Npoints )
    z = zeros( Float64, Npoints )
    #
    nabla2_x = Lmat*x
    r = b - nabla2_x
    z = apply_prec_ilu0( prec, r)
    p = copy(z)

    rsold = dot( r, z )

    for iter = 1 : NiterMax
        #
        nabla2_x = Lmat*p
        #
        alpha = rsold/dot( p, nabla2_x )
        #
        x = x + alpha * p
        r = r - alpha * nabla2_x
        z = apply_prec_ilu0( prec, r)
        #
        rsnew = dot( z, r )
        # deltars = rsold - rsnew
        if showprogress
            @printf("%8d %18.10e\n", iter, sqrt(abs(rsnew)))
        end
        #
        if sqrt(abs(rsnew)) < TOL
            if convmsg
                @printf("#Convergence achieved in linsolve_pcg: %8d iterations.\n", iter)
            end
            break
        end
        #
        p = z + (rsnew/rsold) * p
        #
        rsold = rsnew
    end
    #
    return x
    #
end # of function

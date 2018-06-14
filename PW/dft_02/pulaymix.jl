function pulaymix!(vin::Array{Float64,1}, vout::Array{Float64,1},
                   beta::Float64, dvmat, vmat, iter::Int64, mixdim::Int64)
    # output [vnew,df,dv]
    Npoints = size(vin)[1]
    #
    # construct and modify vmat and dvmat
    #
    if iter <= mixdim
        vmat[:,iter] = vin
        dvmat[:,iter] = vout - vin
    else
        # delete the first leading column  ???
        #vmat[:,1] = []; dvmat[:,1] = []
        # append vin and vout-vin at the end of vmat and dvmat;
        vmat[:,mixdim] = vin
        dvmat[:,mixdim] = vout - vin
    end

    thresh = 1e-6
    if iter > 1
        ibeg = 2
        iend = min(iter,mixdim)
        A = dvmat[:,ibeg:iend]
        b = dvmat[:,ibeg-1]
        ncols = size(A)[2]
        A = b*ones(1,ncols) - A
        U,S,V = svd(A,thin=true)
        while S[ncols] <  thresh*S[1]
            ibeg = ibeg + 1
            if ibeg > iend 
                @printf("warning: no mixing...\n")
                vnew = vout
                return
            end
            A = dvmat[:,ibeg:iend]
            b = dvmat[:,ibeg-1]
            ncols = size(A)[2]
            A = b*ones(1,ncols) - A
            U,S,V = svd(A,thin=true)
        end
        g2 = V*(Matrix(Diagonal(S))\(U'*b))
        g1 = 1 - sum(g2)
        g = [g1;g2]
        vnew = ( vmat[:,ibeg-1:iend] + beta*dvmat[:,ibeg-1:iend] )*g
    else
        vnew = vout
    end
    return vnew
end

const ⊗ = kron

"""
Build Laplacian matrix using Kronecker products
"""
function get_Laplacian3d_kron( LF::LF3dGrid )
    IIx = speye(LF.Nx)
    IIy = speye(LF.Ny)
    IIz = speye(LF.Nz)
    return LF.LFx.D2jl ⊗ IIy ⊗ IIz +
           IIx ⊗ LF.LFy.D2jl ⊗ IIz +
           IIx ⊗ IIy ⊗ LF.LFz.D2jl
end


"""
A naive implementation to build Laplacian by manually, by searching for non-zero
indices and their values.
"""
function get_Laplacian3d_naive( LF::LF3dGrid )
    N = LF.Nx * LF.Ny * LF.Nz

    NNZ = calc_NNZ(LF)
    Tnz = [false, false, false]
    # Nabla2 matrix
    values = zeros(Float64,NNZ)
    column = zeros(Int64,NNZ)
    row    = zeros(Int64,NNZ)
    inz    = 0
    for ip2 = 1:N^3
        for ip1 = ip2:N^3
            Tnz[:] = false
            i1 = LF.lin2xyz[1,ip1]
            i2 = LF.lin2xyz[1,ip2]
            #
            j1 = LF.lin2xyz[2,ip1]
            j2 = LF.lin2xyz[2,ip2]
            #
            k1 = LF.lin2xyz[3,ip1]
            k2 = LF.lin2xyz[3,ip2]
            #
            nabla2 = 0.0
            #
            if j1 == j2 && k1 == k2
                nabla2 = nabla2 + LF.LFx.D2jl[i1,i2]
                Tnz[1] = true
            end
            #
            if i1 == i2 && k1 == k2
                nabla2 = nabla2 + LF.LFy.D2jl[j1,j2]
                Tnz[2] = true
            end
            #
            if i1 == i2 && j1 == j2
                nabla2 = nabla2 + LF.LFz.D2jl[k1,k2]
                Tnz[3] = true
            end
            #
            if any(Tnz) # update non-zero element
                inz = inz + 1
                # Notice that we do not exploit symmetry here
                values[inz] = nabla2
                column[inz] = ip1
                row[inz]    = ip2
                #@printf("%5d %5d: Non zero found, inz = %8d\n", ip1, ip2, )
                if ip1 != ip2 # not a diagonal element
                    inz = inz + 1 # update NNZ again
                    values[inz] = nabla2
                    column[inz] = ip2
                    row[inz]    = ip1
                end
            end # if Tnz
        end # for ip1
    end # for ip2
    return sparse(row, column, values)
end

function calc_NNZ( LF::LF3dGrid ; verbose=false )
    Npoints = LF.Nx * LF.Ny * LF.Nz
    #
    NNZ = 0
    Tnz = [false, false, false]
    # Nabla2 matrix
    for ip2 = 1:Npoints
        for ip1 = ip2:Npoints
            Tnz[:] = false
            i1 = LF.lin2xyz[1,ip1]
            i2 = LF.lin2xyz[1,ip2]
            #
            j1 = LF.lin2xyz[2,ip1]
            j2 = LF.lin2xyz[2,ip2]
            #
            k1 = LF.lin2xyz[3,ip1]
            k2 = LF.lin2xyz[3,ip2]
            #
            nabla2 = 0.0
            #
            if (j1 == j2 && k1 == k2) || (i1 == i2 && k1 == k2) || (i1 == i2 && j1 == j2)
                NNZ = NNZ + 1
                #@printf("%5d %5d: Non zero found, inz = %8d\n", ip1, ip2, NNZ)
                if ip1 != ip2 # not a diagonal element
                    NNZ = NNZ + 1
                    #@printf("%5d %5d: Non zero found, inz = %8d\n", ip1, ip2, NNZ)
                end
            end
        end # ip1
    end # ip2
    if verbose
        @printf("NNZ = %8d, nonzero percentage = %7.2f%%\n", NNZ, NNZ/(Npoints^2) * 100.0 )
    end
    return NNZ
end


function old_calc_NNZ( LF::LF3dGrid ; verbose=false )
    Npoints = LF.Nx * LF.Ny * LF.Nz
    #
    NNZ = 0
    Tnz = [false, false, false]
    # Nabla2 matrix
    for ip2 = 1:Npoints
        for ip1 = ip2:Npoints
            Tnz[:] = false
            i1 = LF.lin2xyz[1,ip1]
            i2 = LF.lin2xyz[1,ip2]
            #
            j1 = LF.lin2xyz[2,ip1]
            j2 = LF.lin2xyz[2,ip2]
            #
            k1 = LF.lin2xyz[3,ip1]
            k2 = LF.lin2xyz[3,ip2]
            #
            nabla2 = 0.0
            #
            if j1 == j2 && k1 == k2
                Tnz[1] = true
            end
            #
            if i1 == i2 && k1 == k2
                Tnz[2] = true
            end
            #
            if i1 == i2 && j1 == j2
                Tnz[3] = true
            end
            if any(Tnz)
                NNZ = NNZ + 1
                #@printf("%5d %5d: Non zero found, inz = %8d\n", ip1, ip2, NNZ)
                if ip1 != ip2 # not a diagonal element
                    NNZ = NNZ + 1
                    #@printf("%5d %5d: Non zero found, inz = %8d\n", ip1, ip2, NNZ)
                end
            end
        end # ip1
    end # ip2
    if verbose
        @printf("NNZ = %8d, nonzero percentage = %7.2f%%\n", NNZ, NNZ/(Npoints^2) * 100.0 )
    end
    return NNZ
end

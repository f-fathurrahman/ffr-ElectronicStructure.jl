function solve_poisson_FFT( Gvec::GvectorsT, rhoR )
    #
    G2 = Gvec.G2
    Ns = Gvec.Ns
    Npoints = prod(Ns)
    #
    ctmp = 4.0*pi*R_to_G( Ns, rhoR )
    #
    for ip = 2:Npoints
        ctmp[ip] = ctmp[ip]/G2[ip]
    end
    ctmp[1] = 0.0
    return real( G_to_R(Ns,ctmp) )
end

# In case we forget to convert the input, we convert it in this version
function R_to_G( Ns::Array{Int,1}, fR_::Array{Float64,1} )
    fR = convert(Array{ComplexF64,1},fR_)
    out = reshape( fft( reshape(fR,Ns[1],Ns[2],Ns[3]) ), size(fR) )
end

function G_to_R( Ns::Array{Int,1}, fG::Array{ComplexF64,1} )
    out = reshape( ifft( reshape(fG,Ns[1],Ns[2],Ns[3]) ), size(fG) )
end

function G_to_R( Ns::Array{Int,1}, fG )
    out = reshape( ifft( reshape(fG,Ns[1],Ns[2],Ns[3]) ), size(fG) )
end


function c_G_to_R(Ns::Array{Int,1}, fG)
    out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_inv_fft3d, "../libs/fft3d.so"), Void,
           (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
            fG, out, Ns[3],Ns[2],Ns[1] )
    return out
end


function R_to_G( Ns::Array{Int,1}, fR )
    #
    out = reshape( fft( reshape(fR,Ns[1],Ns[2],Ns[3]) ),
          size(fR) )
    #
end


function c_R_to_G( Ns::Array{Int,1}, fR::Array{Complex128,1} )
    #
    out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_fw_fft3d, "../libs/fft3d.so"), Void,
           (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
            fR, out, Ns[3],Ns[2],Ns[1] )
    return out
    #
end


#
# Input array is of type Float64. This input is converted
# to Complex128 before supplied to `fftw_fw_fft3d`
#
function c_R_to_G( Ns::Array{Int,1}, fR::Array{Float64,1} )
    out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
    ccall( (:fftw_fw_fft3d, "../libs/fft3d.so"), Void,
           (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64),
            convert(Array{Complex128,1},fR), out, Ns[3],Ns[2],Ns[1] )
    return out
end

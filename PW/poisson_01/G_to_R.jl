function G_to_R( Ns::Array{Int,1}, psi )
  out = reshape( ifft( reshape(psi,Ns[1],Ns[2],Ns[3]) ),
        size(psi) )
end

function c_G_to_R(Ns::Array{Int,1}, fG)
  out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
  ccall( (:fftw_inv_fft3d, "../libs/fft3d.so"), Void,
    (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64), fG, out,
    Ns[3],Ns[2],Ns[1] )
  return out
end

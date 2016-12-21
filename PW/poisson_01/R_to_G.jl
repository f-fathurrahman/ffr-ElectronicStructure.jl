function R_to_G( Ns::Array{Int,1}, psi )
  out = reshape( fft( reshape(psi,Ns[1],Ns[2],Ns[3]) ),
        size(psi) )
end

function c_R_to_G(Ns::Array{Int,1}, fR)
  out = zeros(Complex128,Ns[1]*Ns[2]*Ns[3])
  ccall( (:fftw_fw_fft3d, "../libs/fft3d.so"), Void,
    (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64), fR, out,
    Ns[3],Ns[2],Ns[1] )
  return out
end

"""
function R_to_G( Ns::Array{Int,1}, psi )
  tmp = reshape(psi,Ns[1],Ns[2],Ns[3])
  plan = plan_fft( tmp )
  out = reshape( plan*tmp, size(psi) )
end
"""

function c_R_to_G(Ns::Array{Int,1}, fR)
  out = zeros(Complex128,Ns[1],Ns[2],Ns[3])
  ccall( (:fftw_fw_fft3d, "fft3d"), Void,
    (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64), fR, out,
    Ns[3],Ns[2],Ns[1] )
  return out
end

function c_G_to_R(Ns::Array{Int,1}, fG)
  out = zeros(Complex128,Ns[1],Ns[2],Ns[3])
  ccall( (:fftw_inv_fft3d, "fft3d"), Void,
    (Ptr{Complex128}, Ptr{Complex128}, Int64, Int64, Int64), fG, out,
    Ns[3],Ns[2],Ns[1] )
  return out
end

# generic, one-column version
function R_to_G( Ns::Array{Int,1}, fR )
  out = reshape( fft( reshape(fR,Ns[1],Ns[2],Ns[3]) ),
                   size(fR) )
end

function G_to_R( Ns::Array{Int,1}, fG )
  out = reshape( ifft( reshape(fG,Ns[1],Ns[2],Ns[3]) ),
                 size(fG) )
  return out
end

function test_main()
  Ns = [150,150,150]
  in1 = rand(Complex128,Ns[1],Ns[2],Ns[3])

  @time out1c = c_R_to_G(Ns,in1)
  @time in2c = c_G_to_R(Ns,out1c)

  @time out1 = R_to_G(Ns,in1)
  @time in2 = G_to_R(Ns,out1)

  println("diff out1 = ", sum( abs(out1c - out1) ))
  println("diff out1 = ", sum( abs(in2c - in2) ))
end

test_main()


function orthonormalize!( LF::LF3dGrid, v )
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  ortho_gram_schmidt!( v )
  v[:,:] = (1/sqrt(ΔV))*v  # the [:,:] is significant
end

function orthonormalize( LF::LF3dGrid, v )
  ΔV = LF.LFx.h * LF.LFy.h * LF.LFz.h
  return ortho_gram_schmidt( v )/sqrt(ΔV)
end

FFTW.set_num_threads(2);
blas_set_num_threads(2);

include("PWGrid.jl")
include("cI.jl")
include("cJ.jl")
include("apply_K.jl")
include("ortho_gram_schmidt.jl")

function test_main()
  #
  const Ns = [30, 30, 30]
  const LatVecs = diagm( [6.0, 6.0, 6.0] )
  #
  pw_grid = PWGrid( Ns, LatVecs )
  const Npoints = pw_grid.Npoints
  const Ω = pw_grid.Ω
  #
  const Nstates = 4
  psi  = randn(Npoints,Nstates) + im*randn(Npoints,Nstates)
  # orthonormalization, according to: dot(psi,psi)*Omega = 1
  psi = scale( 1.0/sqrt(Ω), psi*inv(sqrtm(psi'*psi)) )
  Kpsi = apply_K( pw_grid, psi )
end

test_main()
test_main()

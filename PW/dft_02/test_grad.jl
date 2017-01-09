include("../common/PWGrid_v02.jl")
include("../common/ortho_gram_schmidt.jl")
include("../common/wrappers_fft.jl")

include("EnergiesT.jl")
include("PotentialsT.jl")
include("gen_dr.jl")
include("init_pot_harm_3d.jl")
include("apply_K.jl")
include("apply_V_loc.jl")
include("apply_H.jl")
include("gradE.jl")
include("calc_rho.jl")
include("solve_poisson.jl")
include("LDA_VWN.jl")

function printMatrix( A::Array{Complex128,2} )
  Nrows = size(A)[1]
  Ncols = size(A)[2]
  for ir = 1:Nrows
    for ic = 1:Ncols
      @printf("(%6.3f,%6.3f) ", real(A[ir,ic]), imag(A[ir,ic]))
    end
    @printf("\n")
  end
end

function printMatrix( A::Array{Float64,2} )
  Nrows = size(A)[1]
  Ncols = size(A)[2]
  for ir = 1:Nrows
    for ic = 1:Ncols
      @printf("%6.3f ", A[ir,ic])
    end
    @printf("\n")
  end
end


function ∇E( pw::PWGrid, Potentials, Focc, W::Array{Complex128,2} )

  Ngwx    = size(W)[1]
  Nstates = size(W)[2]
  Ω = pw.Ω
  Ns = pw.Ns
  #
  grad = zeros( Complex128, Ngwx, Nstates )

  F = diagm(Focc)
  HW = apply_H( pw, Potentials, W )

  U = W' * W
  U_sqrt = sqrtm( inv(U) )

  # BbbH
  ℍ = U_sqrt * W' * HW * U_sqrt

  #printMatrix(U)
  println("\nU_sqrt = "); printMatrix(U_sqrt)

  HFH = ℍ*F - F*ℍ

  println("\nℍ ="); printMatrix(ℍ)
  println("\nℍ*F - F*ℍ ="); printMatrix(HFH)

  # Calculation of Q
  mu, V = eig(U)
  denom = sqrt(mu)*ones(1,length(mu))
  denom = denom + denom'

  ℚ = V * ( ( V' * HFH * V ) ./ denom ) * V'

  println("mu = "); println(mu)
  println("V = "); printMatrix(V)
  println("denom = "); printMatrix(denom)
  println("ℚ = "); printMatrix(ℚ)

  grad = (HW - W * W'*HW)*F + W*ℚ

  return grad

end



function test_main( Ns )

  const LatVecs = 6.0*diagm( ones(3) )

  pw = PWGrid( Ns, LatVecs )

  const Ω  = pw.Ω
  const r  = pw.r
  const G  = pw.gvec.G
  const G2 = pw.gvec.G2
  const Npoints = prod(Ns)
  const Ngwx = pw.gvecw.Ngwx

  @printf("Ns   = (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
  @printf("Ngwx = %d\n", Ngwx)

  const actual = Npoints/Ngwx
  const theor = 1/(4*pi*0.25^3/3)
  @printf("Compression: actual, theor: %f , %f\n", actual, theor)

  #
  # Generate array of distances
  #
  center = sum(LatVecs,2)/2
  dr = gen_dr( r, center )
  #
  # Setup potential
  #
  V_ionic = init_pot_harm_3d( pw, dr )
  print("sum(Vpot)*Ω/Npoints = $(sum(V_ionic)*Ω/Npoints)\n");
  #
  const Nstates = 4
  #Focc = 2.0*ones(Nstates)
  Focc = [2.0, 1, 1, 1]

  srand(1234)
  psi = randn(Ngwx,Nstates) + im*randn(Ngwx,Nstates)
  psi = ortho_gram_schmidt(psi)

  Potentials = PotentialsT( V_ionic, zeros(Npoints), zeros(Npoints) )
  rho = calc_rho( pw, Focc, psi )
  Potentials.Hartree = real( G_to_R( Ns, solve_poisson(pw, rho) ) )
  Potentials.XC = excVWN( rho ) + rho .* excpVWN( rho )

  grad2 = gradE( pw, Potentials, Focc, psi )

  grad1 = ∇E( pw, Potentials, Focc, psi )

  println("sum(grad1) = ", sum(grad1))
  println("sum(grad2) = ", sum(grad2))

end

test_main( [30, 30, 30] )

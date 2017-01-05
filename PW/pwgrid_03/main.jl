include("../common/PWGrid_v02.jl")

function test_main()
  Ns = [32, 32, 32]
  LatVecs = 8.0*diagm(ones(3))

  pw = PWGrid( Ns, LatVecs )
  @printf("Sampling points: (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])

  actual = prod(pw.Ns)/pw.gvecw.Ngwx
  theor = 1/(4*pi*0.25^3/3)
  @printf("Actual, theor: %10.5f %10.5f\n", actual, theor)

end

@time test_main()

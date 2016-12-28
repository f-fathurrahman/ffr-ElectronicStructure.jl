include("../common/PWGrid_v02.jl")

function test_main()
  const ecutwfc = 35 * 0.5  # in Ha
  LatVecs = 8.0*diagm(ones(3))

  pw = PWGrid( ecutwfc, LatVecs )

  actual = prod(pw.Ns)/pw.gvecw.Ngwx
  theor = 1/(4*pi*0.25^3/3)
  @printf("Actual, theor: %10.5f %10.5f\n", actual, theor)

end

@time test_main()

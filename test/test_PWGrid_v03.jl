using ElectronicStructure

function test_main()
  LatVecs = 10.0*diagm(ones(3))
  pw = PWGrid( 20.0, LatVecs)
  Ns = pw.Ns
  @printf("Sampling points: (%d,%d,%d)\n", Ns[1], Ns[2], Ns[3])
  actual = prod(pw.Ns)/pw.gvecw.Ngwx
  theor = 1/(4*pi*0.25^3/3)
  @printf("Actual, theor: %10.5f %10.5f\n", actual, theor)
end

test_main()

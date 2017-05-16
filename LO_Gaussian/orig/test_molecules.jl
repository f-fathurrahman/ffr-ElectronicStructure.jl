function test_h2()
  @time Energy, E, U = rhf(h2)
  Eref = -1.1170996
  #@assert isapprox(Energy,-1.1170996)
  @printf("Error = %18.10f\n", abs(Energy - Eref))
end

function test_lih()
  @time Energy, E, U = rhf(lih)
  Eref = -7.86073270525799
  #@assert isapprox(Energy,-7.86073270525799)
  @printf("Error = %18.10f\n", abs(Energy - Eref))
end

function test_h2o()
  @time Energy,E,U = rhf(h2o)
  #@assert isapprox(Energy,-74.9597609118851)
  Eref = -74.9597609118851
  @printf("Error = %18.10f\n", abs(Energy - Eref))
end

function test_ch4()
  @time Energy,E,U = rhf(ch4)
  @printf("Energy = %18.10f\n", Energy)
end

function test_c6h6()
  @time Energy,E,U = rhf(c6h6)
  @printf("Energy = %18.10f\n", Energy)
end
